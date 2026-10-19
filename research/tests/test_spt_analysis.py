"""Verify quantitative SPT reductions and the complete resumable project workflow."""

import json
import xml.etree.ElementTree as ET
from dataclasses import replace
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pytest
import xarray as xr

from research.spt.analysis import main, process_analysis
from research.spt.analysis_data import (
    AnalysisOptions,
    CellSampler,
    load_analysis,
    spatial_profile,
    timelines,
    weighted_statistics,
)
from research.spt.analysis_plots import PlotOptions
from research.spt.data_processing import process_studies
from research.spt.geometry import CONVERGENCE_MESHES, process_geometries
from research.spt.pipeline import main as pipeline_main
from research.spt.pipeline import run_pipeline
from research.spt.plot_data import plot_sources
from research.tests.test_spt_geometry import _source


def _input(
    root: Path,
    name: str,
    count: int,
    times: tuple[float, ...],
    *,
    changing: bool = False,
) -> None:
    path = _source(root, name, count, times)
    tree = ET.parse(path)
    with h5py.File(path.with_suffix(".h5"), "r+") as handle:
        for grid in tree.findall(".//Grid[@CollectionType='Temporal']/Grid"):
            time = grid.find("Time")
            assert time is not None
            fraction = float(time.attrib["Value"]) / times[-1]
            for attribute in grid.findall("Attribute"):
                item = attribute.find("DataItem")
                assert item is not None and item.text
                key = attribute.attrib["Name"]
                if key == "rr_position":
                    handle[item.text.partition(":")[2]][...] = np.linspace(
                        0, 1, count
                    ).reshape(-1, 1)
                elif key == "rr_protein":
                    handle[item.text.partition(":")[2]][...] = np.linspace(
                        0.2, 0.8, count
                    ).reshape(-1, 1)
                elif changing and key in {"rr_(T)", "element_volume_point_TPM"}:
                    values = (
                        [10 * fraction, 10 * (1 - fraction)]
                        if key == "rr_(T)"
                        else [1 + 2 * fraction, 3 - 2 * fraction]
                    )
                    handle[item.text.partition(":")[2]][...] = np.array(values).reshape(
                        -1, 1
                    )


@pytest.fixture
def study(tmp_path: Path) -> tuple[Path, Path]:
    """Provide changing cell volumes that expose summary-interpolation errors."""
    root, output = tmp_path / "raw", tmp_path / "local"
    _input(root, "sim001", 2, (0, 4), changing=True)
    assert (
        process_geometries(
            root,
            output,
            study="simulations",
            run_ids=("sim001",),
            combinations=False,
            workers=1,
        )[0].status
        == "processed"
    )
    assert (
        process_studies(
            root,
            output,
            study="simulations",
            run_ids=("sim001",),
            geometry="lobulus",
            workers=1,
        )[0].status
        == "processed"
    )
    return root, output


def test_cell_first_statistics_and_serialized_results(
    study: tuple[Path, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    """Interpolate cells and weights before reducing; never reopen network HDF data."""
    root, output = study
    monkeypatch.setattr(
        "h5py.File",
        lambda *args, **kwargs: pytest.fail("Analysis must use local numerical caches"),
    )
    result = process_analysis(
        root,
        output,
        output,
        study="simulations",
        run_ids=("sim001",),
        options=AnalysisOptions(time_points=3, bins=4),
        plots=PlotOptions(formats=("svg", "png"), dpi=45),
        workers=1,
    )
    assert [item.status for item in result] == ["processed"]
    source = plot_sources(
        root, output, study="simulations", geometry="lobulus", run_ids=("sim001",)
    )[0]
    with (
        load_analysis(source, output, "native") as native,
        load_analysis(source, output, "aligned") as aligned,
        load_analysis(source, output, "profiles") as profile,
        load_analysis(source, output, "cells") as cells,
    ):
        np.testing.assert_allclose(native["rr_(T)_mean"], [7.5, 7.5])
        assert float(aligned["mean"].sel(time=2, field="rr_(T)")) == 5
        assert float(aligned["mean"].sel(time=2, field="rr_necrosis")) == 0
        assert float(aligned["mean"].sel(time=4, field="rr_necrosis")) == 1
        assert np.isnan(profile["mean"].isel(position=1)).all()
        np.testing.assert_array_equal(profile.cell_count.isel(time=0), [6, 0, 0, 6])
        np.testing.assert_allclose(cells["rr_(T)"].sel(time=2), 5)
        assert list(cells.time.values) == [0, 2, 4]
    directory = output / "analysis/simulations_lobulus"
    manifest = json.loads((directory / "manifest.json").read_text())
    generation = directory / manifest["generation"]
    with xr.open_dataset(generation / "analysis.nc", engine="scipy") as data:
        assert list(data.run.values) == ["sim001"]
        assert (
            float(data.aligned_mean.sel(run="sim001", time=2, analysis_field="rr_(T)"))
            == 5
        )
    native_table = pd.read_csv(generation / "native-timecourses.tsv", sep="\t")
    assert list(native_table.time) == [0, 4]
    assert "Necrotic volume [%]" in (generation / "timecourses.svg").read_text()
    assert not list(
        ET.parse(generation / "timecourses.svg").iter(
            "{http://www.w3.org/2000/svg}image"
        )
    )
    assert len(list(generation.glob("*.png"))) == 6


def test_unchanged_analysis_is_metadata_only(
    study: tuple[Path, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    """Reuse data and figures without numerical reads or worker startup."""
    root, output = study
    options, plots = (
        AnalysisOptions(time_points=3, bins=4),
        PlotOptions(formats=("svg",)),
    )
    process_analysis(
        root,
        output,
        output,
        study="simulations",
        run_ids=("sim001",),
        options=options,
        plots=plots,
        workers=1,
    )

    def unexpected(*args: object, **kwargs: object) -> None:
        pytest.fail("Unchanged analysis must not open numerical data or start workers")

    for name in (
        "xarray.open_dataset",
        "research.spt.analysis.ProcessPoolExecutor",
        "research.spt.analysis_plots.load_analysis",
    ):
        monkeypatch.setattr(name, unexpected)
    result = process_analysis(
        root,
        output,
        output,
        study="simulations",
        run_ids=("sim001",),
        options=options,
        plots=plots,
        workers=2,
    )
    assert result[0].status == "cached"
    report = json.loads((output / "analysis-report-simulations.json").read_text())
    assert report["native_frames_read"] == 0


def test_empty_bins_endpoints_and_invalid_weights() -> None:
    """Keep endpoint cells, preserve fractional necrosis, and reject invalid physics."""
    frame = {
        "rr_position": np.array([0.0, 1.0]),
        "element_volume_point_TPM": np.array([1.0, 3.0]),
        "rr_necrosis": np.array([0.25, 1.0]),
    }
    profile = spatial_profile(frame, ("rr_necrosis",), 4)
    np.testing.assert_array_equal(profile.cell_count, [1, 0, 0, 1])
    assert float(profile["mean"].isel(field=0, position=0)) == 0.25
    assert np.isnan(profile["mean"].isel(field=0, position=1))
    with pytest.raises(ValueError, match="volumes"):
        weighted_statistics(np.array([[1, 2]]), np.array([0, 0]))
    with pytest.raises(ValueError, match="rr_position"):
        spatial_profile(
            {**frame, "rr_position": np.array([-0.2, 1.0])}, ("rr_necrosis",), 4
        )
    with pytest.raises(ValueError, match="finite"):
        weighted_statistics(np.array([[np.nan, 2]]), np.array([1, 1]))
    dataset = xr.Dataset(
        {"rr_necrosis": (("time", "cell"), [[0.25, 0], [1, 1]])},
        coords={"time": [0, 2], "cell": [0, 1]},
    )
    np.testing.assert_array_equal(
        CellSampler(dataset, ("rr_necrosis",)).sample(1)["rr_necrosis"], [0.25, 0]
    )


def test_failure_recovery_and_missing_artifact(
    study: tuple[Path, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    """Keep previous figures when a replacement fails and repair missing exports."""
    root, output = study
    options, plots = (
        AnalysisOptions(time_points=3, bins=4),
        PlotOptions(formats=("svg",)),
    )
    process_analysis(
        root,
        output,
        output,
        study="simulations",
        run_ids=("sim001",),
        options=options,
        plots=plots,
        workers=1,
    )
    directory = output / "analysis/simulations_lobulus"
    manifest = directory / "manifest.json"
    previous = manifest.read_bytes()
    generation = directory / json.loads(previous)["generation"]
    with monkeypatch.context() as patch:

        def fail(*args: object, **kwargs: object) -> None:
            raise RuntimeError("Synthetic rendering failure")

        patch.setattr("research.spt.analysis_plots._curves", fail)
        results = process_analysis(
            root,
            output,
            output,
            study="simulations",
            run_ids=("sim001",),
            options=options,
            plots=plots,
            workers=1,
            force=True,
        )
        assert any(item.status == "failed" for item in results)
    assert manifest.read_bytes() == previous
    assert (generation / "timecourses.svg").is_file()
    (generation / "timecourses.svg").unlink()
    results = process_analysis(
        root,
        output,
        output,
        study="simulations",
        run_ids=("sim001",),
        options=options,
        plots=plots,
        workers=1,
    )
    assert all(item.status == "cached" for item in results)
    current = directory / json.loads(manifest.read_text())["generation"]
    assert (current / "timecourses.svg").is_file()
    assert not list(output.rglob(".partial-*"))


def test_full_pipeline_from_fresh_inputs(tmp_path: Path) -> None:
    """Run every stage from raw sixths through mixed views, analysis and figures."""
    root, output = tmp_path / "raw", tmp_path / "local"
    for index, key in enumerate((*CONVERGENCE_MESHES, "50000")):
        _input(root, f"convergence_{key}", 8 - index, (0, 2) if index == 0 else (0, 4))
    before = {
        p: (p.stat().st_size, p.stat().st_mtime_ns)
        for p in root.rglob("*")
        if p.is_file()
    }
    report = run_pipeline(
        root,
        output,
        study="convergence",
        workers=2,
        formats=("svg", "png"),
        analysis_options=AnalysisOptions(time_points=3, bins=4),
        plot_options=PlotOptions(
            formats=("svg", "png"), dpi=45, convergence_symlog=False
        ),
        geometry_snapshots=1,
    )
    assert report["status"] == "passed", report
    assert [item["name"] for item in report["stages"]] == [
        "catalogue",
        "geometry",
        "numeric-cache",
        "analysis",
        "geometry-figures",
    ]
    assert (output / "index.html").is_file()
    assert (output / "analysis/index.html").is_file()
    assert (
        len(
            list((output / "analysis-data/convergence_lobulus").glob("*/manifest.json"))
        )
        == 7
    )
    directory = output / "analysis/convergence_lobulus"
    generation = (
        directory / json.loads((directory / "manifest.json").read_text())["generation"]
    )
    differences = pd.read_csv(generation / "convergence-means.tsv", sep="\t")
    assert set(differences.reference_run) == {"convergence_06250"}
    assert differences.time.max() == 2
    assert (
        differences.absolute_difference[differences.run == "convergence_06250"]
        .eq(0)
        .all()
    )
    native = pd.read_csv(generation / "native-timecourses.tsv", sep="\t")
    assert native[native.run == "convergence_50000"].time.max() == 4
    assert {p: (p.stat().st_size, p.stat().st_mtime_ns) for p in before} == before
    assert (
        pipeline_main(
            [
                "--input-root",
                str(root),
                "--output-root",
                str(output),
                "--study",
                "convergence",
                "--status",
                "--time-points",
                "3",
                "--bins",
                "4",
                "--geometry-snapshots",
                "1",
                "--formats",
                "svg",
                "png",
            ]
        )
        == 0
    )


def test_dry_run_and_invalid_selections(
    study: tuple[Path, Path], tmp_path: Path
) -> None:
    """Keep inspection read-only and reject unsupported analysis inputs."""
    root, output = study
    destination = tmp_path / "new-analysis"
    args = [
        "--input-root",
        str(root),
        "--data-root",
        str(output),
        "--output-root",
        str(destination),
        "--study",
        "simulations",
        "--runs",
        "sim001",
    ]
    assert main([*args, "--status"]) == 0
    assert not destination.exists()
    assert main([*args, "--bins", "0"]) == 1
    assert main([*args, "--spatial-times", "5"]) == 1
    assert main([*args, "--fields", "pressure"]) == 1
    assert main([*args, "--runs", "dose_01"]) == 1
    report = run_pipeline(
        root, tmp_path / "empty", study="simulations", run_ids=("sim001",), dry_run=True
    )
    assert report["status"] == "pending"
    assert not (tmp_path / "empty").exists()


def test_pipeline_stops_after_input_failure(tmp_path: Path) -> None:
    """Avoid downstream writes and record the failed stage clearly."""
    output = tmp_path / "out"
    report = run_pipeline(tmp_path / "missing", output, study="convergence")
    assert report["status"] == "failed"
    assert len(report["stages"]) == 1
    assert not (output / "analysis").exists()
    assert (
        json.loads((output / "pipeline-report.json").read_text())["status"] == "failed"
    )


def test_later_spatial_times_keep_initial_diagnostics(
    study: tuple[Path, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    """Keep the actual initial frame and reuse reductions when only styling changes."""
    root, output = study
    options = AnalysisOptions(
        fields=("rr_(T)",), time_points=3, bins=4, spatial_times=(2, 4)
    )
    result = process_analysis(
        root,
        output,
        output,
        study="simulations",
        run_ids=("sim001",),
        options=options,
        plots=PlotOptions(formats=("svg",), dpi=45),
        workers=1,
    )
    assert result[0].status == "processed"
    source = plot_sources(
        root, output, study="simulations", geometry="lobulus", run_ids=("sim001",)
    )[0]
    with (
        load_analysis(source, output, "initial") as initial,
        load_analysis(source, output, "profiles") as profiles,
    ):
        assert list(initial.time.values) == [0]
        assert list(profiles.time.values) == [2, 4]
        assert float(initial["mean"].sel(field="rr_(T)").isel(time=0, position=0)) == 0

    def unexpected(*args: object, **kwargs: object) -> None:
        pytest.fail("Display changes must reuse numerical reductions")

    monkeypatch.setattr("research.spt.analysis_data.CellSampler.sample", unexpected)
    result = process_analysis(
        root,
        output,
        output,
        study="simulations",
        run_ids=("sim001",),
        options=options,
        plots=PlotOptions(formats=("svg",), dpi=60),
        workers=1,
    )
    assert result[0].status == "cached"


@pytest.mark.filterwarnings("error:overflow encountered:RuntimeWarning")
def test_zero_difference_symlog_is_finite() -> None:
    """Render zero differences without inventing negative errors or overflowing."""
    from research.spt.analysis_plots import _scale
    from research.spt.plot_render import figure

    art = figure((4, 3))
    try:
        axis = art.subplots()
        axis.plot([1, 2], [0, 0])
        _scale(axis, "difference", True)
        art.canvas.draw()
        assert axis.get_ylim()[0] == 0
        assert np.isfinite(axis.get_yticks()).all()
    finally:
        art.clear()


def test_status_lists_pending_stages_after_geometry(
    study: tuple[Path, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Inspect a partially built project without HDF reads or new files."""
    root, _ = study
    output = tmp_path / "geometry-only"
    process_geometries(
        root,
        output,
        study="simulations",
        run_ids=("sim001",),
        combinations=False,
        workers=1,
    )
    before = set(output.rglob("*"))
    monkeypatch.setattr(
        "h5py.File",
        lambda *args, **kwargs: pytest.fail("Status must not open HDF data"),
    )
    report = run_pipeline(
        root, output, study="simulations", run_ids=("sim001",), dry_run=True
    )
    assert report["status"] == "pending"
    assert [item["name"] for item in report["stages"]] == [
        "catalogue",
        "geometry",
        "numeric-cache",
        "analysis",
        "geometry-figures",
    ]
    assert set(output.rglob("*")) == before


def test_study_reference_times_are_preserved(study: tuple[Path, Path]) -> None:
    """Retain 10 h simulation and 1.5 h convergence comparisons within valid bounds."""
    root, output = study
    source = plot_sources(
        root, output, study="simulations", geometry="lobulus", run_ids=("sim001",)
    )[0]
    simulation = replace(source, end=86400)
    convergence = replace(
        source, run=replace(source.run, study="convergence_lobulus"), end=36590
    )
    windows = timelines([simulation, convergence], AnalysisOptions(time_points=3))
    assert windows["simulations"][1] == (0, 36000, 86400)
    assert windows["convergence"][1] == (0, 5400, 36590)
