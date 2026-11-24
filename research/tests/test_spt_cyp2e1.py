"""Verify the CYP2E1 necrosis exports and their figures."""

import json
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from research.spt.analysis import process_analysis
from research.spt.analysis_cyp2e1 import cyp2e1_frames, render_cyp2e1_figures
from research.spt.analysis_data import AnalysisOptions
from research.spt.analysis_plots import PlotOptions
from research.spt.catalog import RunSource
from research.spt.conditions import SimulationCondition, boundary_flows
from research.spt.data_processing import process_studies
from research.spt.geometry import process_geometries
from research.spt.plot_data import PlotSource
from research.spt.plot_render import RenderOptions
from research.tests.test_spt_analysis import _input

FINAL_TIME = 86400.0
PROFILE_FIELDS = ("rr_necrosis", "rr_protein")


def _profiles(
    patterns: list[int],
    doses: list[int],
    bins: int,
    *,
    empty_bin: int | None = None,
) -> tuple[xr.Dataset, list[PlotSource]]:
    """Build a profile dataset and its sources without touching any solver output."""
    times = (0.0, FINAL_TIME)
    edges = np.linspace(0.0, 1.0, bins + 1)
    centers = (edges[:-1] + edges[1:]) / 2
    runs, sources = [], []
    for index, (pattern, dose) in enumerate(zip(patterns, doses, strict=True), start=1):
        run_id = f"sim{index:03}"
        runs.append(run_id)
        condition = SimulationCondition(run_id, pattern, dose, boundary_flows[dose])
        source = RunSource(
            Path("raw"), "simulations_lobulus", run_id, condition, "test"
        )
        sources.append(
            PlotSource(source, "lobulus", 0.0, FINAL_TIME, run_id, dose=dose)
        )
    shape = (len(runs), len(times), len(PROFILE_FIELDS), bins)
    mean = np.tile(np.linspace(0.1, 0.9, bins), (*shape[:-1], 1))
    count = np.full((len(runs), len(times), bins), 5, dtype=np.int32)
    if empty_bin is not None:
        mean[..., empty_bin] = np.nan
        count[..., empty_bin] = 0
    data = {
        "mean": (("run", "time", "field", "position"), mean),
        "std": (("run", "time", "field", "position"), mean / 10),
        "minimum": (("run", "time", "field", "position"), mean / 2),
        "maximum": (("run", "time", "field", "position"), mean * 2),
        "cell_count": (("run", "time", "position"), count),
        "volume": (("run", "time", "position"), count.astype(float) * 1e-12),
        "bin_left": (
            ("run", "time", "position"),
            np.tile(edges[:-1], (len(runs), len(times), 1)),
        ),
        "bin_right": (
            ("run", "time", "position"),
            np.tile(edges[1:], (len(runs), len(times), 1)),
        ),
    }
    profiles = xr.Dataset(
        data,
        coords={
            "run": list(runs),
            "time": list(times),
            "field": list(PROFILE_FIELDS),
            "position": centers,
        },
    )
    return profiles, sources


def test_no_cyp2e1_runs_produce_no_frames() -> None:
    """Skip the export for a selection that holds no CYP2E1 run."""
    profiles, sources = _profiles(patterns=[0], doses=[0], bins=4)
    necrosis, paired = cyp2e1_frames(profiles, sources)
    assert necrosis.empty and paired.empty


def test_necrosis_frame_columns_and_rows() -> None:
    """Name the flux of every row and keep one row per position bin."""
    profiles, sources = _profiles(patterns=[6, 6], doses=[0, 1], bins=4)
    necrosis, _ = cyp2e1_frames(profiles, sources)
    assert list(necrosis.columns) == [
        "run",
        "pattern",
        "pattern_name",
        "dose",
        "boundary_flow",
        "time",
        "position",
        "bin_left",
        "bin_right",
        "necrosis_mean",
        "necrosis_std",
        "necrosis_minimum",
        "necrosis_maximum",
        "cell_count",
        "volume",
    ]
    assert len(necrosis) == 8
    assert set(necrosis.pattern_name) == {"cyp2e1"}
    assert list(necrosis.dose.unique()) == [1, 2]
    assert necrosis.time.unique().tolist() == [FINAL_TIME]


def test_source_without_condition_is_ignored() -> None:
    """Drop a run whose catalogue row carries no condition instead of failing."""
    profiles, sources = _profiles(patterns=[6, 6], doses=[0, 1], bins=4)
    sources[1] = replace(sources[1], run=replace(sources[1].run, condition=None))
    necrosis, paired = cyp2e1_frames(profiles, sources)
    assert list(necrosis.run.unique()) == [sources[0].run.run_id]
    assert list(paired.run.unique()) == [sources[0].run.run_id]


def test_paired_frame_pairs_protein_and_necrosis_in_one_bin() -> None:
    """Pair both quantities of one position bin in a single row."""
    profiles, sources = _profiles(patterns=[6], doses=[0], bins=4)
    _, paired = cyp2e1_frames(profiles, sources)
    assert list(paired.columns) == [
        "run",
        "pattern",
        "pattern_name",
        "dose",
        "boundary_flow",
        "time",
        "position",
        "bin_left",
        "bin_right",
        "protein_mean",
        "protein_std",
        "necrosis_mean",
        "necrosis_std",
    ]
    assert len(paired) == 4


def test_empty_bins_keep_their_rows() -> None:
    """Keep a bin without cells, so position indexing is never silently shortened."""
    profiles, sources = _profiles(patterns=[6], doses=[0], bins=4, empty_bin=1)
    necrosis, paired = cyp2e1_frames(profiles, sources)
    assert len(necrosis) == 4 and len(paired) == 4
    assert int(necrosis.cell_count.iloc[1]) == 0
    assert np.isnan(necrosis.necrosis_mean.iloc[1])
    assert np.isnan(paired.protein_mean.iloc[1])


def test_figures_are_written_in_every_format(tmp_path: Path) -> None:
    """Write both figures in each requested format."""
    profiles, sources = _profiles(patterns=[6, 6], doses=[0, 1], bins=4)
    necrosis, paired = cyp2e1_frames(profiles, sources)
    written = render_cyp2e1_figures(
        necrosis, paired, tmp_path, RenderOptions(formats=("svg", "png"), dpi=45)
    )
    assert sorted(path.name for path in written) == [
        "cyp2e1-necrosis-position.png",
        "cyp2e1-necrosis-position.svg",
        "cyp2e1-necrosis-protein.png",
        "cyp2e1-necrosis-protein.svg",
    ]
    assert all(path.is_file() for path in written)


def test_figure_labels_and_legend(tmp_path: Path) -> None:
    """Name both axes and every substrate flux in the legend."""
    profiles, sources = _profiles(patterns=[6, 6], doses=[0, 7], bins=4)
    necrosis, paired = cyp2e1_frames(profiles, sources)
    render_cyp2e1_figures(
        necrosis, paired, tmp_path, RenderOptions(formats=("svg",), dpi=45)
    )
    document = (tmp_path / "cyp2e1-necrosis-position.svg").read_text()
    assert "D1" in document and "D8" in document
    assert "Necrosis" in document and "position" in document.lower()
    assert "CYP2E1" in (tmp_path / "cyp2e1-necrosis-protein.svg").read_text()


def test_single_run_and_empty_frames(tmp_path: Path) -> None:
    """Render a one-flux selection and write nothing at all for empty frames."""
    profiles, sources = _profiles(patterns=[6], doses=[3], bins=4)
    necrosis, paired = cyp2e1_frames(profiles, sources)
    written = render_cyp2e1_figures(
        necrosis, paired, tmp_path, RenderOptions(formats=("png",), dpi=45)
    )
    assert len(written) == 2
    assert (
        render_cyp2e1_figures(
            pd.DataFrame(),
            pd.DataFrame(),
            tmp_path,
            RenderOptions(formats=("png",), dpi=45),
        )
        == []
    )


def _cyp2e1_study(tmp_path: Path, run_ids: tuple[str, ...]) -> tuple[Path, Path]:
    """Prepare CYP2E1 runs through geometry and numerical caches."""
    root, output = tmp_path / "raw", tmp_path / "local"
    for run_id in run_ids:
        _input(root, run_id, 2, (0, 4), changing=True)
    assert all(
        item.status == "processed"
        for item in process_geometries(
            root,
            output,
            study="simulations",
            run_ids=run_ids,
            combinations=False,
            workers=1,
        )
    )
    assert all(
        item.status == "processed"
        for item in process_studies(
            root,
            output,
            study="simulations",
            run_ids=run_ids,
            geometry="lobulus",
            workers=1,
        )
    )
    return root, output


def test_analysis_publishes_cyp2e1_exports(tmp_path: Path) -> None:
    """Publish both tables and both figures inside the analysis generation."""
    root, output = _cyp2e1_study(tmp_path, ("sim049", "sim050"))
    assert [
        item.status
        for item in process_analysis(
            root,
            output,
            output,
            study="simulations",
            run_ids=("sim049", "sim050"),
            options=AnalysisOptions(time_points=3, bins=4),
            plots=PlotOptions(formats=("svg", "png"), dpi=45),
            workers=1,
        )
    ] == ["processed", "processed"]
    directory = output / "analysis/simulations_lobulus"
    generation = (
        directory / json.loads((directory / "manifest.json").read_text())["generation"]
    )
    table = pd.read_csv(generation / "cyp2e1-necrosis-position.tsv", sep="\t")
    assert len(table) == 8
    assert sorted(table.dose.unique()) == [1, 2]
    assert set(table.pattern_name) == {"cyp2e1"}
    paired = pd.read_csv(generation / "cyp2e1-necrosis-protein.tsv", sep="\t")
    assert {"protein_mean", "necrosis_mean"} <= set(paired.columns)
    for name in ("cyp2e1-necrosis-position", "cyp2e1-necrosis-protein"):
        assert (generation / f"{name}.svg").is_file()
        assert (generation / f"{name}.png").is_file()


def test_selection_without_cyp2e1_writes_no_exports(tmp_path: Path) -> None:
    """Write no CYP2E1 artifact at all for a selection without that pattern."""
    root, output = _cyp2e1_study(tmp_path, ("sim001",))
    process_analysis(
        root,
        output,
        output,
        study="simulations",
        run_ids=("sim001",),
        options=AnalysisOptions(time_points=3, bins=4),
        plots=PlotOptions(formats=("png",), dpi=45),
        workers=1,
    )
    directory = output / "analysis/simulations_lobulus"
    generation = (
        directory / json.loads((directory / "manifest.json").read_text())["generation"]
    )
    assert not list(generation.glob("cyp2e1-*"))


def test_analysis_catalogue_records_the_export_paths(tmp_path: Path) -> None:
    """Name both published tables in the stage catalogue."""
    root, output = _cyp2e1_study(tmp_path, ("sim049", "sim050"))
    process_analysis(
        root,
        output,
        output,
        study="simulations",
        run_ids=("sim049", "sim050"),
        options=AnalysisOptions(time_points=3, bins=4),
        plots=PlotOptions(formats=("png",), dpi=45),
        workers=1,
    )
    catalogue = pd.read_csv(output / "logs/analysis-simulations.tsv", sep="\t")
    assert catalogue.cyp2e1_position.iloc[0].endswith("cyp2e1-necrosis-position.tsv")
    assert catalogue.cyp2e1_protein.iloc[0].endswith("cyp2e1-necrosis-protein.tsv")
