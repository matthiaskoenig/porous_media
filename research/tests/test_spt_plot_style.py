"""Protect the manuscript's scientific display choices against styling drift."""

from itertools import pairwise
from pathlib import Path

import numpy as np
import pytest
import xarray as xr
from matplotlib.collections import PolyCollection
from matplotlib.container import ErrorbarContainer
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle

from research.spt.analysis_plots import PlotOptions, _curves, _diagnostics, _resolution
from research.spt.catalog import RunSource
from research.spt.conditions import simulation_conditions
from research.spt.fields import LAYERS, POSITION, SPT_LAYERS, VOLUME
from research.spt.plot_data import PlotSource, sample_times
from research.spt.plot_geometry import Prepared, shared_ranges
from research.spt.plot_render import RenderOptions, Samples, render_timecourse
from research.spt.plot_style import convergence_color, simulation_color
from research.spt.visualization import data_layers_spt


def _source(root: Path, *, convergence: bool = False) -> PlotSource:
    """Describe a plot without requiring simulation files or private data."""
    run = RunSource(
        root,
        "convergence_lobulus" if convergence else "simulations_lobulus",
        "convergence_06250" if convergence else "sim016",
        None if convergence else simulation_conditions()[15],
        "synthetic",
    )
    return PlotSource(
        run, "lobulus", 0, 86400, run.run_id, dose=None if convergence else 7
    )


def _cells() -> xr.Dataset:
    """Keep unequal cells, fractional necrosis and declared diagnostic units."""
    return xr.Dataset(
        {
            POSITION: (("time", "cell"), [[0, 0.5, 1]] * 2),
            "rr_protein": (("time", "cell"), [[1, 2, 3]] * 2),
            "rr_(T)": (("time", "cell"), [[0, 2, 4], [1, 3, 5]]),
            "rr_necrosis": (("time", "cell"), [[0.25, 0.5, 1]] * 2),
            VOLUME: (("time", "cell"), [[1e-12, 2e-12, 3e-12]] * 2),
            "rr_Vext": (("time", "cell"), [[1e-9, 2e-9, 3e-9]] * 2),
            "rr_Vli": (("time", "cell"), [[2e-9, 3e-9, 4e-9]] * 2),
            "pressure": (("time", "cell"), [[1, 2, 3]] * 2),
        },
        coords={"time": [0, 5400], "cell": [0, 1, 2]},
    )


def test_canonical_mesh_layers_and_legacy_adapter() -> None:
    """Keep the original labels and colormaps in both rendering entry points."""
    assert [(layer.sid, layer.colormap) for layer in LAYERS] == [
        ("rr_(S_ext)", "magma"),
        ("rr_(P_ext)", "magma"),
        ("rr_protein", "hot"),
        ("rr_(T)", "magma"),
        ("rr_necrosis", "binary"),
    ]
    assert [layer.title for layer in LAYERS] == [
        "Substrate S plasma [mM]",
        "Product P plasma [mM]",
        "Protein",
        "Toxic compound T [mM]",
        "Necrosis (0: alive, 1: death)",
    ]
    for canonical, legacy in zip(SPT_LAYERS, data_layers_spt, strict=True):
        assert canonical.sid == legacy.sid
        assert canonical.colormap == legacy.colormap
        assert canonical.title == legacy.title
        assert canonical.limits == legacy.color_limits is None
        assert canonical.scalar_bar == legacy.scalar_bar is True


def test_original_geometry_snapshots_and_full_history_limits(tmp_path: Path) -> None:
    """Use the five manuscript rows and include extrema between sampled frames."""
    source = _source(tmp_path)
    chosen = sample_times([source], 5, manuscript_snapshots=True)["simulations"][0]
    np.testing.assert_allclose(chosen, np.linspace(0, 86400, 10)[1::2])
    assert sample_times([source], 5)["simulations"][0][0] == 0
    assert sample_times([source], 5, (0, 3600), manuscript_snapshots=True)[
        "simulations"
    ][0] == (0, 3600)
    descriptor = source.run.xdmf
    descriptor.parent.mkdir(parents=True, exist_ok=True)
    descriptor.with_name(f"{descriptor.stem}_limits.json").write_text(
        '{"limits": {"rr_(T)": [0, 100]}}'
    )
    prepared = Prepared(source, {"metadata": {"ranges": {"rr_(T)": [2, 5]}}})
    assert shared_ranges([prepared])["simulations"]["rr_(T)"] == (0, 100)


def test_local_limits_override_malformed_legacy_metadata(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Read validated local metadata instead of failing on empty raw sidecars."""
    source = _source(tmp_path)
    source.run.xdmf.parent.mkdir(parents=True)
    source.run.xdmf.with_name(f"{source.run.xdmf.stem}_limits.json").write_text("")
    cache = tmp_path / "cache"
    cache.mkdir()
    (cache / "limits.json").write_text(
        '{"Cell": {"rr_(T)": {"min": 0, "max": 42, "nonfinite": 0}}}'
    )
    monkeypatch.setattr(
        "research.spt.plot_geometry._current_generation", lambda *args: cache
    )
    prepared = Prepared(source, {"metadata": {"ranges": {"rr_(T)": [2, 5]}}})
    assert shared_ranges([prepared], tmp_path)["simulations"]["rr_(T)"] == (0, 42)
    assert shared_ranges([prepared])["simulations"]["rr_(T)"] == (2, 5)
    source.run.xdmf.with_name(f"{source.run.xdmf.stem}_limits.json").write_text(
        '{"limits": {"rr_(T)": [0, 100], "rr_(S)": [0, 8]}}'
    )
    prepared = Prepared(
        source, {"metadata": {"ranges": {"rr_(T)": [2, 5], "rr_(S)": [1, 2]}}}
    )
    ranges = shared_ranges([prepared], tmp_path)["simulations"]
    assert ranges["rr_(T)"] == (0, 42)
    assert ranges["rr_(S)"] == (0, 8)


def test_mesh_edges_and_per_panel_scalar_bars(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Retain cell boundaries and a labeled five-tick horizontal bar in every panel."""
    data = Samples(
        np.array([[[0, 0], [1, 0], [0, 1]]]),
        np.array([0]),
        np.array([[[0.5]], [[1.0]]]),
        np.array([-1]),
        np.array([0, 1]),
    )
    layer = next(layer for layer in LAYERS if layer.sid == "rr_protein")

    def inspect(art: Figure, destination: Path, options: RenderOptions) -> None:
        assert tuple(art.get_size_inches()) == (4, 8)
        assert len(art.axes) == 4
        for axis in art.axes[:2]:
            collection = axis.collections[0]
            assert isinstance(collection, PolyCollection)
            assert collection.cmap.name == "hot"
            np.testing.assert_allclose(
                np.asarray(collection.get_edgecolor(), dtype=float)[0],
                [169 / 255] * 3 + [1],
            )
            assert np.asarray(collection.get_linewidth(), dtype=float)[0] == 0.48
        for axis in art.axes[2:]:
            assert axis.get_title() == "Protein"
            assert len(axis.get_xticks()) == 5

    monkeypatch.setattr("research.spt.plot_render.save_figure", inspect)
    render_timecourse(
        _source(tmp_path),
        data,
        tmp_path,
        RenderOptions(layers=(layer,), snapshots=(0, 1)),
        {layer.sid: (0, 1)},
    )


def test_small_colorbar_labels_remain_distinct() -> None:
    """Do not turn nonzero product concentrations into five identical zero labels."""
    from research.spt.plot_render import figure, scalar_bar

    art = figure((4, 4))
    try:
        axis = art.subplots()
        collection = PolyCollection([[[0, 0], [1, 0], [0, 1]]])
        collection.set_clim(0, 0.006)
        scalar_bar(art, axis, collection, LAYERS[1])
        art.set_dpi(150)
        art.canvas.draw()
        labels = [tick.get_text() for tick in art.axes[-1].get_xticklabels()]
        assert len(set(labels)) == 5
        exponent = art.axes[-1].xaxis.get_offset_text().get_text().replace("−", "-")
        assert "10^{-3}" in exponent
        boxes = [tick.get_window_extent() for tick in art.axes[-1].get_xticklabels()]
        assert not any(left.overlaps(right) for left, right in pairwise(boxes))
    finally:
        art.clear()


def test_simulation_errorbars_and_cell_scatter(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Keep dose colors, SD error bars, bold labels and every fractional cell value."""
    source = _source(tmp_path)
    aligned = xr.Dataset(
        {
            "mean": (("time", "field"), [[1], [2]]),
            "std": (("time", "field"), [[0.1], [0.2]]),
        },
        coords={"time": [0, 86400], "field": ["rr_(T)"]},
    )

    def inspect(art: Figure, destination: Path, options: RenderOptions) -> None:
        axis = art.axes[0]
        assert tuple(art.get_size_inches()) == (2.7, 2.5)
        assert axis.yaxis.label.get_fontweight() == "bold"
        assert axis.yaxis.label.get_fontsize() == 12
        assert not art.legends
        assert not axis.get_title()
        assert axis.lines[0].get_color() == "#08306bff"
        assert axis.lines[0].get_marker() == "o"
        assert axis.lines[0].get_markeredgecolor() == "black"
        if destination.name == "timecourses":
            assert isinstance(axis.containers[0], ErrorbarContainer)
            assert axis.containers[0].has_yerr
            assert axis.get_xlabel() == "Time [hr]"
            np.testing.assert_allclose(axis.get_ylim(), [-0.11, 2.31])
        else:
            assert axis.get_ylabel() == "Necrosis [%]"
            np.testing.assert_allclose(
                np.asarray(axis.lines[0].get_ydata(), dtype=float), [25, 50, 100]
            )
            assert [tick.get_text() for tick in axis.get_xticklabels()] == [
                "PP",
                "",
                "PV",
            ]

    monkeypatch.setattr("research.spt.analysis_plots.save_figure", inspect)
    _curves([source], [aligned], tmp_path, "timecourses", PlotOptions(), ("rr_(T)",))
    _curves(
        [source],
        [_cells()],
        tmp_path,
        "spatial",
        PlotOptions(),
        ("rr_necrosis",),
        spatial_time=5400,
    )
    assert simulation_color(1, 8) == "#08306bff"
    assert convergence_color("06250", 0, 1) == "#fde725"
    assert convergence_color("50000", 0, 1) == "#440154"


def test_initial_histogram_and_diagnostic_conversions(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Show cell counts and original nL conversions instead of binned summaries."""

    def inspect(art: Figure, destination: Path, options: RenderOptions) -> None:
        assert all(isinstance(patch, Rectangle) for patch in art.axes[0].patches)
        assert (
            sum(
                patch.get_bbox().height
                for patch in art.axes[0].patches
                if isinstance(patch, Rectangle)
            )
            == 3
        )
        assert len(art.axes[0].patches) == 10
        assert art.axes[0].get_ylabel() == "Frequency [-]"
        assert art.axes[2].get_ylabel() == "Mesh cell volume [nl]"
        np.testing.assert_allclose(
            np.asarray(art.axes[2].lines[0].get_ydata(), dtype=float), [1, 2, 3]
        )
        np.testing.assert_allclose(
            np.asarray(art.axes[3].lines[0].get_ydata(), dtype=float), [1, 2, 3]
        )
        np.testing.assert_allclose(
            np.asarray(art.axes[4].lines[0].get_ydata(), dtype=float), [2, 3, 4]
        )

    monkeypatch.setattr("research.spt.analysis_plots.save_figure", inspect)
    _diagnostics([_source(tmp_path)], [_cells().isel(time=0)], tmp_path, PlotOptions())


def test_convergence_lines_and_resolution_reference(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Keep the reference line, upper legend and linear resolution axis at 1.5 h."""
    source = _source(tmp_path, convergence=True)
    aligned = xr.Dataset(
        {
            "mean": (("time", "field"), [[1], [2]]),
            "std": (("time", "field"), [[0.1], [0.2]]),
        },
        coords={"time": [0, 86400], "field": ["rr_(T)"]},
    )

    def inspect(art: Figure, destination: Path, options: RenderOptions) -> None:
        axis = art.axes[0]
        if destination.name == "timecourses":
            assert axis.lines[0].get_marker() == "None"
            assert isinstance(axis.containers[0], ErrorbarContainer)
            assert not axis.containers[0].has_yerr
            np.testing.assert_allclose(
                np.asarray(axis.lines[1].get_xdata(), dtype=float), [1.5, 1.5]
            )
            assert axis.lines[1].get_linestyle() == "--"
            assert len(art.legends) == 1
            np.testing.assert_allclose(axis.get_ylim(), [0, 2.05])
        else:
            assert axis.get_xscale() == "linear"
            assert axis.get_xlabel() == "# Elements"
            # Protein is the first field in this synthetic selection.
            np.testing.assert_allclose(
                np.asarray(art.axes[1].lines[0].get_ydata(), dtype=float), [22 / 6]
            )

    monkeypatch.setattr("research.spt.analysis_plots.save_figure", inspect)
    _curves([source], [aligned], tmp_path, "timecourses", PlotOptions(), ("rr_(T)",))
    _resolution(
        [_cells()], [source], tmp_path, PlotOptions(), False, ("rr_protein", "rr_(T)")
    )
