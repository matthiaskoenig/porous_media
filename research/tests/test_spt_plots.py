"""Exercise selective SPT figures, cache recovery and optional streamed videos."""

import json
import shutil
import subprocess
import xml.etree.ElementTree as ET
from dataclasses import replace
from pathlib import Path
from typing import cast

import h5py
import meshio
import numpy as np
import pytest

from porous_media.visualization.pyvista_visualization import DataLayer
from research.spt.cache import ProcessingResult
from research.spt.catalog import discover_runs
from research.spt.geometry import CONVERGENCE_MESHES, process_geometries
from research.spt.plot_data import (
    PlotSource,
    plot_sources,
    prepare_samples,
    projected_faces,
    sample_times,
)
from research.spt.plot_geometry import main, process_plots
from research.spt.plot_render import Layer, RenderOptions, Samples, display_limits
from research.spt.plots_simulations_geometry import visualize_spt_2d
from research.tests.test_spt_geometry import _source

LAYERS = (
    Layer("rr_(T)", "Toxicant"),
    Layer("rr_necrosis", "Necrosis", "Reds", (0, 1)),
)
OPTIONS = RenderOptions(LAYERS, ("svg", "pdf", "png"), dpi=45)


@pytest.fixture
def simulation(tmp_path: Path) -> tuple[Path, Path]:
    """Supply a small extruded mesh with continuous and discrete fields."""
    root, output = tmp_path / "raw", tmp_path / "output"
    _source(root, "sim001", 2, (0.0, 2.0, 4.0))
    return root, output


def _samples(source: PlotSource, output: Path) -> tuple[Samples, dict[str, object]]:
    manifest = prepare_samples(
        source, output, tuple(layer.sid for layer in LAYERS), (0, 1, 2, 3, 4)
    )
    assert manifest is not None
    path = source.directory(output, "plot-data") / str(manifest["generation"])
    return Samples.load(path / "samples.npz"), cast(
        dict[str, object], manifest["metadata"]
    )


def test_selected_fields_times_and_cell_values(simulation: tuple[Path, Path]) -> None:
    """Ignore unreadable unused arrays and hold necrosis between native outputs."""
    root, output = simulation
    descriptor = root / "simulations_sixth/xdmf/sim001.xdmf"
    tree = ET.parse(descriptor)
    for item in tree.findall(".//Attribute[@Name='flux']/DataItem"):
        item.text = "sim001.h5:/does-not-exist"
    tree.write(descriptor)
    source = plot_sources(root, output, study="simulations", geometry="sixth")[0]
    data, metadata = _samples(source, output)
    assert metadata["hdf_field_reads"] == 6
    assert data.values.shape == (5, 2, 2)
    np.testing.assert_array_equal(data.frame(1)[0], [4, 4])
    np.testing.assert_array_equal(data.frame(1)[1], [0, 0])
    np.testing.assert_array_equal(data.frame(3)[1], [1, 1])
    np.testing.assert_array_equal(data.cell_ids, [0, 1])
    assert data.faces.shape == (2, 3, 2)
    with pytest.raises(ValueError, match="outside source history"):
        prepare_samples(source, output, ("rr_(T)",), (5,))
    with pytest.raises(ValueError, match="Missing scalar"):
        prepare_samples(source, output, ("flux",), (0,))


def test_planar_quads_and_invalid_projection() -> None:
    """Triangulate a quad without averaging cells and reject tilted extrusions."""
    points = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]], dtype=float)
    faces, cell_ids = projected_faces(
        points, [meshio.CellBlock("quad", [[0, 1, 2, 3]])]
    )
    assert faces.shape == (2, 3, 2)
    np.testing.assert_array_equal(cell_ids, [0, 0])
    wedge = np.vstack((points[:3], points[:3] + np.array([0.1, 0, 1])))
    with pytest.raises(ValueError, match=r"vertical|extrusion|planar"):
        projected_faces(wedge, [meshio.CellBlock("wedge", [[0, 1, 2, 3, 4, 5]])])
    with pytest.raises(ValueError, match="Unsupported"):
        projected_faces(points, [meshio.CellBlock("tetra", [[0, 1, 2, 3]])])


def test_vector_figures_and_unchanged_resume(
    simulation: tuple[Path, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    """Resume without XML, HDF, mesh reads, sample loading or worker startup."""
    root, output = simulation
    before = {
        path: (path.stat().st_size, path.stat().st_mtime_ns)
        for path in root.rglob("*")
        if path.is_file()
    }
    results = process_plots(
        root,
        output,
        output,
        study="simulations",
        geometry="sixth",
        options=OPTIONS,
        snapshots=2,
        workers=1,
    )
    assert [result.status for result in results] == ["rendered"]
    row = results[0].row
    svg = ET.parse(str(row["svg"]))
    assert not list(svg.iter("{http://www.w3.org/2000/svg}image"))
    assert list(svg.iter("{http://www.w3.org/2000/svg}path"))
    assert Path(str(row["pdf"])).read_bytes().startswith(b"%PDF")
    assert Path(str(row["png"])).read_bytes().startswith(b"\x89PNG")
    assert not list(output.rglob("*.xdmf"))
    assert not list(output.rglob("*.mp4"))
    assert not list(output.rglob("panels"))
    assert "SVG" in (output / "figures/index-simulations.html").read_text()

    def unexpected(*args: object, **kwargs: object) -> None:
        pytest.fail("Unchanged plots must use metadata only")

    for name in (
        "research.spt.plot_data.inspect_run",
        "meshio.xdmf.TimeSeriesReader",
        "h5py.File",
        "research.spt.plot_render.Samples.load",
        "research.spt.execution.ProcessPoolExecutor",
    ):
        monkeypatch.setattr(name, unexpected)
    resumed = process_plots(
        root,
        output,
        output,
        study="simulations",
        geometry="sixth",
        options=OPTIONS,
        snapshots=2,
        workers=2,
    )
    assert [result.status for result in resumed] == ["cached"]
    report = json.loads((output / "plot-report-simulations.json").read_text())
    assert report["sampling"]["hdf_field_reads"] == 0
    assert {
        path: (path.stat().st_size, path.stat().st_mtime_ns) for path in before
    } == before


def test_repair_and_failed_publication(
    simulation: tuple[Path, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    """Repair a removed artifact and keep the previous generation on render failure."""
    root, output = simulation

    def run(force: bool = False) -> list[ProcessingResult]:
        return process_plots(
            root,
            output,
            output,
            study="simulations",
            geometry="sixth",
            options=OPTIONS,
            snapshots=1,
            workers=1,
            overviews=False,
            force=force,
        )

    first = run()[0]
    png = Path(str(first.row["png"]))
    png.unlink()
    repaired = run()[0]
    assert repaired.status == "rendered"
    assert Path(str(repaired.row["png"])).is_file()
    manifest = output / "figures/simulations_sixth/sim001/manifest.json"
    before = manifest.read_bytes()

    def fail(*args: object) -> None:
        raise RuntimeError("Synthetic save failure")

    monkeypatch.setattr("research.spt.plot_render.render_timecourse", fail)
    assert run(force=True)[0].status == "failed"
    assert manifest.read_bytes() == before
    assert Path(str(repaired.row["svg"])).is_file()
    assert not list(output.rglob(".partial-*"))


def test_three_geometries_and_separate_common_windows(tmp_path: Path) -> None:
    """Render all available convergence meshes and the six-resolution comparison."""
    root, output = tmp_path / "raw", tmp_path / "local"
    for index, key in enumerate((*CONVERGENCE_MESHES, "50000")):
        _source(
            root,
            f"convergence_{key}",
            index + 1,
            (0, 1, 2) if index == 0 else (0, 2, 4),
        )
    built = process_geometries(root, output, study="convergence", workers=1)
    assert all(result.status == "processed" for result in built)
    sources = plot_sources(root, output, study="convergence")
    assert len(sources) == 15
    mixed = next(source for source in sources if source.geometry == "combined")
    assert mixed.sector_labels == CONVERGENCE_MESHES
    assert sample_times(sources, 1)["convergence"][0] == (2,)
    simulation = replace(
        sources[-1],
        run=replace(sources[-1].run, study="simulations_sixth"),
        start=0,
        end=24,
    )
    windows = sample_times([*sources, simulation], 2)
    assert windows["simulations"][0] == (0, 24)
    assert windows["convergence"][0] == (0, 2)
    results = process_plots(
        root,
        output,
        output,
        study="convergence",
        options=replace(OPTIONS, formats=("svg", "png")),
        snapshots=2,
        workers=1,
    )
    assert len(results) == 15
    assert all(result.status == "rendered" for result in results)
    catalog = (output / "plot-catalog-convergence.tsv").read_text()
    assert "convergence_combinations" in catalog
    assert len(list((output / "overviews").rglob("*.svg"))) == 6


def test_generator_compatibility_and_selected_panels(
    simulation: tuple[Path, Path],
) -> None:
    """Retain custom layer scales and accept generators without rendering other fields."""
    root, output = simulation
    paths = (source.xdmf for source in discover_runs(root, ("simulations_sixth",)))
    layer = DataLayer("rr_(T)", "Custom units", conversion_factor=-2)
    results = visualize_spt_2d(
        paths,
        [layer, DataLayer("not_selected", "Ignore")],
        output,
        (name for name in ["rr_(T)"]),
        create_panels=True,
        snapshots=1,
        formats=("svg",),
    )
    assert results[0].status == "rendered"
    panels = list(output.rglob("panels/*.svg"))
    assert len(panels) == 1
    assert "Custom units" in panels[0].read_text()
    assert layer.color_limits is None
    assert display_limits(Layer("T", "T", factor=-2), {"T": (1, 4)}) == (-8, -2)
    with pytest.raises(ValueError, match="at least one"):
        visualize_spt_2d(iter(()), [layer], output, [layer.sid])


def test_dry_run_invalid_options_and_overviews_only(
    simulation: tuple[Path, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    """Validate selections before writing; ordinary figures do not require FFmpeg."""
    root, output = simulation
    args = [
        "--input-root",
        str(root),
        "--data-root",
        str(output),
        "--output-root",
        str(output),
        "--study",
        "simulations",
        "--geometry",
        "sixth",
    ]
    assert main([*args, "--status"]) == 0
    assert not output.exists()
    assert main([*args, "--workers", "0"]) == 1
    assert main([*args, "--runs", "missing"]) == 1
    assert main([*args, "--times", "5"]) == 1
    assert main([*args, "--only-overviews", "--no-overviews"]) == 1
    monkeypatch.setattr("research.spt.plot_render.shutil.which", lambda name: None)
    results = process_plots(
        root,
        output,
        output,
        study="simulations",
        geometry="sixth",
        options=OPTIONS,
        snapshots=1,
        workers=1,
        per_history=False,
    )
    assert results[0].status == "sampled"
    assert not list((output / "figures").rglob("timecourse.*"))
    assert len(list((output / "overviews").rglob("*.svg"))) == 2
    with pytest.raises(FileNotFoundError, match="FFmpeg"):
        process_plots(
            root,
            output,
            output,
            study="simulations",
            geometry="sixth",
            options=OPTIONS,
            video_frames=3,
            workers=1,
        )


@pytest.mark.skipif(shutil.which("ffmpeg") is None, reason="FFmpeg is unavailable")
def test_streamed_video_has_requested_frames(simulation: tuple[Path, Path]) -> None:
    """Stream three frames with odd DPI and no individual panel or frame images."""
    root, output = simulation
    results = process_plots(
        root,
        output,
        output,
        study="simulations",
        geometry="sixth",
        options=replace(OPTIONS, layers=LAYERS[:1], formats=("svg",), dpi=51, fps=2),
        snapshots=1,
        video_frames=3,
        workers=1,
        overviews=False,
    )
    assert results[0].status == "rendered", results[0].message
    video = Path(str(results[0].row["video"]))
    completed = subprocess.run(
        ["ffmpeg", "-i", str(video), "-f", "null", "-"],
        capture_output=True,
        check=True,
    )
    assert b"frame=    3" in completed.stderr
    assert not list(output.rglob("*.png"))
    assert not list(output.rglob("panels"))


def test_changed_data_invalidate_samples(simulation: tuple[Path, Path]) -> None:
    """Rebuild compact samples after a source change instead of reusing stale values."""
    root, output = simulation
    source = plot_sources(root, output, study="simulations", geometry="sixth")[0]
    original, _ = _samples(source, output)
    with h5py.File(source.run.xdmf.with_suffix(".h5"), "r+") as handle:
        tree = ET.parse(source.run.xdmf)
        item = tree.find(".//Attribute[@Name='rr_(T)']/DataItem")
        assert item is not None and item.text
        handle[item.text.partition(":")[2]][...] = 99
    rebuilt, _ = _samples(source, output)
    assert original.frame(0)[0, 0] == 2
    assert rebuilt.frame(0)[0, 0] == 99


def test_output_paths_cannot_write_into_sources(simulation: tuple[Path, Path]) -> None:
    """Reject overlapping roots and a redirected plot directory before output writes."""
    root, output = simulation
    with pytest.raises(ValueError, match="nonoverlapping"):
        process_plots(
            root,
            output,
            root / "plots",
            study="simulations",
            geometry="sixth",
            workers=1,
        )
    output.mkdir()
    (output / "plot-data").symlink_to(root, target_is_directory=True)
    with pytest.raises(ValueError, match="symbolic links"):
        process_plots(
            root, output, output, study="simulations", geometry="sixth", workers=1
        )
    assert not (root / "simulations_sixth/sim001").exists()


def test_render_options_reuse_samples(
    simulation: tuple[Path, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    """Change display resolution without repeating source array reads."""
    root, output = simulation
    first = process_plots(
        root,
        output,
        output,
        study="simulations",
        geometry="sixth",
        options=OPTIONS,
        snapshots=1,
        workers=1,
        overviews=False,
    )[0]

    def unexpected(*args: object, **kwargs: object) -> None:
        pytest.fail("Rendering options must reuse existing plot samples")

    monkeypatch.setattr("h5py.File", unexpected)
    second = process_plots(
        root,
        output,
        output,
        study="simulations",
        geometry="sixth",
        options=replace(OPTIONS, dpi=65),
        snapshots=1,
        workers=1,
        overviews=False,
    )[0]
    assert second.status == "rendered", second.message
    assert second.row["png"] != first.row["png"]


def test_explicit_symlog_figures(simulation: tuple[Path, Path]) -> None:
    """Label nonlinear scales and preserve the linear necrosis range and source values."""
    root, output = simulation
    assert (
        main(
            [
                "--input-root",
                str(root),
                "--data-root",
                str(output),
                "--output-root",
                str(output),
                "--study",
                "simulations",
                "--geometry",
                "sixth",
                "--scale",
                "symlog",
                "--fields",
                "rr_(T)",
                "rr_necrosis",
                "--snapshots",
                "1",
                "--formats",
                "svg",
                "--no-overviews",
                "--workers",
                "1",
            ]
        )
        == 0
    )
    directory = output / "figures/simulations_sixth/sim001"
    manifest = json.loads((directory / "manifest.json").read_text())
    layers = manifest["configuration"]["options"]["layers"]
    assert layers[0]["scale"] == "symlog"
    assert layers[1]["scale"] == "linear"
    assert layers[1]["limits"] == [0, 1]
    svg = directory / manifest["generation"] / "timecourse.svg"
    assert "Toxicant [mM] (symlog)" in svg.read_text()
    assert not list(ET.parse(svg).iter("{http://www.w3.org/2000/svg}image"))
    report = json.loads((output / "plot-report-simulations.json").read_text())
    assert report["ranges"]["simulations"]["rr_(T)"] == [10, 10]
