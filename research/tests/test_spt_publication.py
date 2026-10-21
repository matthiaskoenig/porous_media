"""Check publication annotations, true sampling times, local reads and cache reuse."""

import json
import xml.etree.ElementTree as ET
from dataclasses import replace
from pathlib import Path

import matplotlib as mpl
import numpy as np
import pytest
import xarray as xr
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
from matplotlib.patches import Patch
from PIL import Image

from research.spt.analysis_data import AnalysisOptions, prepare_analysis
from research.spt.conditions import simulation_conditions
from research.spt.data_processing import process_studies
from research.spt.geometry import CONVERGENCE_MESHES, process_geometries
from research.spt.plot_data import PlotSource, plot_sources, prepare_samples
from research.spt.plot_render import RenderOptions, Samples
from research.spt.publication import (
    create_publication,
    remove_obsolete_compact,
    simulation_cohort,
)
from research.spt.publication_data import (
    CONVERGENCE_DISPLAY_MESHES,
    SIMULATION_TIMES,
    PublicationData,
    convergence_view,
    pattern_label,
)
from research.spt.publication_geometry import (
    raster_panel,
    render_example,
    render_necrosis,
    render_patterns,
    save_publication,
)
from research.spt.publication_patterns import (
    PatternIcon,
    add_pattern_icon,
    pattern_profile,
)
from research.spt.publication_timecourses import (
    CONVERGENCE_COLORS,
    annotate_timecourses,
    render_convergence,
)
from research.spt.report import build_report
from research.tests.test_spt_analysis import _input
from research.tests.test_spt_report import ReportParser, _publish

pytestmark = pytest.mark.filterwarnings("error:constrained_layout not applied")


@pytest.fixture
def publication_data(tmp_path: Path) -> PublicationData:
    """Provide one real synthetic lobule with unequal, changing cell volumes."""
    raw, root = tmp_path / "raw", tmp_path / "results"
    _input(raw, "sim014", 2, (0, 86400), changing=True)
    process_geometries(
        raw,
        root,
        study="simulations",
        run_ids=("sim014",),
        combinations=False,
        workers=1,
    )
    process_studies(
        raw,
        root,
        study="simulations",
        run_ids=("sim014",),
        geometry="lobulus",
        workers=1,
    )
    source = plot_sources(
        raw, root, study="simulations", geometry="lobulus", run_ids=("sim014",)
    )[0]
    assert prepare_analysis(
        source,
        root,
        root,
        AnalysisOptions(time_points=3),
        (0, 43200, 86400),
        (0, 86400),
    )
    assert prepare_samples(source, root, AnalysisOptions().fields, (0, 86400))
    return PublicationData(root)


def _no_hdf(*args: object, **kwargs: object) -> None:
    pytest.fail("Publication figures must not open raw or derived HDF histories")


def test_true_times_and_discrete_necrosis_use_local_caches(
    publication_data: PublicationData,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A 6 h label represents 6 h data, with discrete states held before native changes."""
    monkeypatch.setattr("h5py.File", _no_hdf)
    monkeypatch.setattr("meshio.xdmf.TimeSeriesReader", _no_hdf)
    data = publication_data
    source = data.source("simulations_lobulus", "sim014")
    samples = data.snapshots(source, SIMULATION_TIMES)
    np.testing.assert_array_equal(samples.times / 3600, [0, 6, 12, 18, 24])
    np.testing.assert_allclose(samples.values[1, 3, :2], [2.5, 7.5])
    np.testing.assert_array_equal(samples.values[:-1, -1], 0)
    np.testing.assert_array_equal(samples.values[-1, -1], 1)
    assert source.run.condition and source.run.condition.boundary_flow_key == 5
    assert samples is data.snapshots(source, SIMULATION_TIMES)
    with pytest.raises(ValueError, match="outside native history"):
        data.snapshots(source, (0, 86401))
    with pytest.raises(ValueError, match="six patterns and eight doses"):
        simulation_cohort(data)


def test_map_headers_time_arrow_and_vector_annotations(
    publication_data: PublicationData,
    tmp_path: Path,
) -> None:
    """Publication maps retain editable SVG labels over dense, bounded raster fields."""
    data = publication_data
    source = data.source("simulations_lobulus", "sim014")
    with mpl.rc_context({"svg.fonttype": "none"}):
        render_example(
            data,
            source,
            data.ranges([source]),
            tmp_path,
            RenderOptions(formats=("svg",), dpi=35),
        )
    tree = ET.parse(tmp_path / "figure-04.svg")
    text = " ".join(tree.getroot().itertext())
    assert all(
        label in text
        for label in [
            "Substrate",
            "Product",
            "Protein",
            "Necrosis",
            "Time",
            "0 h",
            "6 h",
            "12 h",
            "18 h",
            "24 h",
        ]
    )
    assert list(tree.iter("{http://www.w3.org/2000/svg}image"))
    assert list(tree.iter("{http://www.w3.org/2000/svg}text"))
    labels = [
        "".join(node.itertext())
        for node in tree.iter("{http://www.w3.org/2000/svg}text")
    ]
    assert labels.count("Substrate S plasma [mM]") == 1
    assert any("B  Linear increase | D6" in label for label in labels)
    assert len(list(tree.iter("{http://www.w3.org/2000/svg}image"))) == 25
    assert (
        len([node for node in tree.iter() if node.get("id") == "pattern-icon-b"]) == 1
    )


def test_map_tiles_do_not_keep_full_composite_raster_buffers(
    publication_data: PublicationData,
) -> None:
    """Each PDF image owns only one tile, even within the 150-panel composite."""
    from research.spt.fields import LAYERS
    from research.spt.plot_render import figure

    data = publication_data
    source = data.source("simulations_lobulus", "sim014")
    samples = data.snapshots(source, (0,))
    art = figure((18, 27))
    try:
        axis = art.subplots()
        mappable = raster_panel(
            axis, samples, LAYERS[0], samples.values[0, 0], (0, 10), side=1.5, dpi=100
        )
        pixels = axis.images[0].get_array()
        assert pixels is not None
        assert 0 < pixels.shape[0] <= 150 and 0 < pixels.shape[1] <= 150
        assert pixels.shape[:2] != (150, 150)
        assert pixels.nbytes < 150 * 150 * 4
        assert not axis.collections
        assert mappable.cmap.name == "magma"
        assert mappable.get_clim() == (0, 10)
    finally:
        art.clear()


def test_pattern_panels_and_necrosis_doses_are_annotated(
    publication_data: PublicationData,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A-F pattern order and the eight-dose necrosis grid follow recorded conditions."""
    data = publication_data
    prototype = data.source("simulations_lobulus", "sim014")
    all_sources = [
        replace(
            prototype,
            run=replace(prototype.run, run_id=c.sim_id, condition=c),
            dose=c.boundary_flow_key,
        )
        for c in simulation_conditions()
    ]
    original = data.snapshots(prototype, SIMULATION_TIMES)

    def snapshots(source: PlotSource, times: tuple[float, ...]) -> Samples:
        return replace(
            original,
            times=np.asarray(times),
            values=original.values[-1:] if len(times) == 1 else original.values,
        )

    monkeypatch.setattr(data, "snapshots", snapshots)
    options = RenderOptions(formats=("svg",), dpi=25)
    with mpl.rc_context({"svg.fonttype": "none"}):
        render_patterns(
            data,
            [s for s in all_sources if s.dose == 5],
            data.ranges([prototype]),
            tmp_path,
            options,
        )
        render_necrosis(data, all_sources, tmp_path, options)
    patterns = (tmp_path / "figure-05.svg").read_text()
    assert all(
        label in patterns for label in [pattern_label(pattern) for pattern in range(6)]
    )
    necrosis = (tmp_path / "figure-07.svg").read_text()
    assert all(f"D{dose}" in necrosis for dose in range(1, 9))
    assert "Dose " not in necrosis
    assert "Substrate / drug inflow" in necrosis
    assert patterns.count("Substrate S plasma [mM]") == 1
    assert "D6" in patterns
    assert all(pattern_label(pattern) in necrosis for pattern in range(6))
    assert necrosis.count(">necrosis<") == 1
    assert necrosis.count(">no necrosis<") == 1
    assert "0 (alive)" not in necrosis and "1 (death)" not in necrosis
    for figure in (patterns, necrosis):
        for pattern in range(6):
            assert figure.count(f'id="pattern-icon-{chr(97 + pattern)}"') == 1


def test_shared_bars_and_tight_geometry_panels(
    publication_data: PublicationData, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Shared scales replace panel bars and map axes reach the geometric tile bounds."""
    data = publication_data
    prototype = data.source("simulations_lobulus", "sim014")
    all_sources = [
        replace(
            prototype,
            run=replace(prototype.run, run_id=c.sim_id, condition=c),
            dose=c.boundary_flow_key,
        )
        for c in simulation_conditions()
    ]
    original = data.snapshots(prototype, SIMULATION_TIMES)

    def snapshots(source: PlotSource, times: tuple[float, ...]) -> Samples:
        return replace(
            original,
            times=np.asarray(times),
            values=original.values[-1:] if len(times) == 1 else original.values,
        )

    monkeypatch.setattr(data, "snapshots", snapshots)

    def inspect(
        art: Figure,
        destination: Path,
        options: RenderOptions,
        *,
        width_mm: float | None = None,
    ) -> None:
        maps = [axis for axis in art.axes if axis.images]
        bars = [axis for axis in art.axes if axis.axison]
        if destination.name == "figure-07":
            assert len(maps) == 48 and not bars
            legend = art.legends[0]
            assert [text.get_text() for text in legend.get_texts()] == [
                "necrosis",
                "no necrosis",
            ]
            assert [
                handle.get_facecolor()
                for handle in legend.legend_handles
                if isinstance(handle, Patch)
            ] == [
                (0, 0, 0, 1),
                (1, 1, 1, 1),
            ]
            art.canvas.draw()
            assert legend.get_window_extent().x0 > max(
                axis.get_window_extent().x1 for axis in maps
            )
        else:
            assert len(maps) == (25 if destination.name == "figure-04" else 150)
            assert len(bars) == 5
        for axis in maps:
            extent = axis.images[0].get_extent()
            np.testing.assert_allclose(axis.get_xlim(), extent[:2])
            np.testing.assert_allclose(axis.get_ylim(), extent[2:])

    monkeypatch.setattr("research.spt.publication_geometry.save_publication", inspect)
    options = RenderOptions(formats=("png",), dpi=25)
    ranges = data.ranges([prototype])
    render_example(data, prototype, ranges, tmp_path, options)
    render_patterns(
        data, [s for s in all_sources if s.dose == 5], ranges, tmp_path, options
    )
    render_necrosis(data, all_sources, tmp_path, options)


def test_timecourse_pattern_and_dose_labels_are_consistent() -> None:
    """Time-course rows use the same pattern names and D1-D8 labels as map figures."""
    from research.spt.plot_render import figure

    art = figure((13.5, 15))
    try:
        axes = art.subplots(6, 5, squeeze=False)
        annotate_timecourses(art, axes)
        expected = [
            "A  Constant",
            "B  Linear increase",
            "C  Sharp perivenous",
            "D  Linear decrease",
            "E  Sharp periportal",
            "F  Random",
        ]
        assert [text.get_text() for text in art.legends[0].get_texts()] == [
            f"D{dose}" for dose in range(1, 9)
        ]
        for row in range(6):
            assert [text.get_text() for text in axes[row, 0].texts] == [expected[row]]
            assert len(axes[row, 0].findobj(match=PatternIcon)) == 1
    finally:
        art.clear()


def test_pattern_icons_show_zonation_and_scale_as_vector_art(tmp_path: Path) -> None:
    """Monochrome profiles distinguish patterns and survive compact vector exports."""
    from research.spt.plot_render import figure

    profiles = []
    for pattern in range(6):
        position, values = pattern_profile(pattern)
        profiles.append((values[0], values[-1]))
        np.testing.assert_array_equal(pattern_profile(pattern)[1], values)
        assert np.isfinite(position).all() and np.isfinite(values).all()
        assert position[0] == 0 and position[-1] == 1
        assert np.all(np.diff(position) >= 0)
    assert profiles[:5] == [(0.5, 0.5), (0, 1), (0, 1), (1, 0), (1, 0)]
    for pattern in (2, 4):
        position, values = pattern_profile(pattern)
        assert position[1] == position[2] and values[1] != values[2]
    assert np.var(pattern_profile(5)[1]) > 0.02
    with pytest.raises(ValueError, match="Unknown zonation pattern"):
        pattern_profile(6)
    art = figure((8, 2))
    try:
        icons = []
        for pattern, axis in enumerate(art.subplots(1, 6)):
            axis.set_axis_off()
            label = axis.text(0.7, 0.5, chr(65 + pattern))
            icon = add_pattern_icon(label, pattern, size=30)
            assert icon.profile.get_color() == "black"
            assert icon.panel.get_facecolor() == (1, 1, 1, 1)
            assert icon.panel.get_edgecolor() == (0, 0, 0, 1)
            icons.append(icon)
        options = RenderOptions(formats=("svg", "png"), dpi=100)
        save_publication(art, tmp_path / "patterns", options)
        canvas = art.canvas
        assert isinstance(canvas, FigureCanvasAgg)
        before = icons[0].get_window_extent(canvas.get_renderer())
        save_publication(art, tmp_path / "patterns-compact", options, width_mm=101.6)
        after = icons[0].get_window_extent(canvas.get_renderer())
        assert after.width == pytest.approx(before.width / 2)
        for stem in ("patterns", "patterns-compact"):
            tree = ET.parse(tmp_path / f"{stem}.svg")
            assert not list(tree.iter("{http://www.w3.org/2000/svg}image"))
            for pattern in range(6):
                assert any(
                    node.get("id") == f"pattern-icon-{chr(97 + pattern)}"
                    for node in tree.iter()
                )
    finally:
        art.clear()


def test_convergence_supplement_preserves_native_ends_and_shows_meshes(
    publication_data: PublicationData,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Mesh-only panels and both axis scales use full, genuinely truncated histories."""
    data = publication_data
    prototype = data.source("simulations_lobulus", "sim014")
    samples = data.mesh(prototype.run.study, prototype.run.run_id)
    sources = [
        replace(
            prototype,
            run=replace(
                prototype.run,
                run_id=f"convergence_{mesh}",
                condition=None,
                study="convergence_lobulus",
            ),
            end={"06250": 36590.0, "09375": 77687.6}.get(mesh, 86400.0),
        )
        for mesh in CONVERGENCE_MESHES
    ]
    monkeypatch.setattr(data, "mesh", lambda study, run_id: samples)
    snapshot = data.snapshots(prototype, (5400.0,))
    monkeypatch.setattr(data, "mixed_snapshot", lambda time: snapshot)
    ranges = data.ranges([prototype])
    monkeypatch.setattr(data, "ranges", lambda sources: ranges)

    def summary(source: PlotSource, *, native: bool = False) -> xr.Dataset:
        assert native
        return xr.Dataset(
            {
                f"{field}_mean": ("time", [0, 0.5, 1])
                for field in AnalysisOptions().fields
            },
            coords={"time": [0, source.end / 2, source.end]},
        )

    monkeypatch.setattr(data, "summary", summary)

    def inspect(
        art: Figure,
        destination: Path,
        options: RenderOptions,
    ) -> None:
        assert destination.name == "figure-s1-convergence"
        assert len(art.axes) == 27
        assert len(art.axes[0].images) == 1
        for column in range(5):
            axis = art.axes[7 + 2 * column]
            assert "at 1.5 h" in axis.get_title()
            assert len(axis.images) == 1
            assert len(axis.collections) == 6
        for row in (17, 22):
            curves = {
                str(line.get_label()): line
                for line in art.axes[row].lines
                if line.get_label() in CONVERGENCE_MESHES
            }
            assert len(curves) == 6
            assert np.asarray(curves["06250"].get_xdata())[-1] == 36590 / 3600
            assert np.asarray(curves["09375"].get_xdata())[-1] == 77687.6 / 3600
            assert np.asarray(curves["37250"].get_xdata())[-1] == 24
            for mesh, curve in curves.items():
                assert curve.get_color() == CONVERGENCE_COLORS[mesh]
        assert [text.get_text() for text in art.legends[0].get_texts()] == list(
            CONVERGENCE_DISPLAY_MESHES
        )
        for index, mesh in enumerate(CONVERGENCE_DISPLAY_MESHES):
            np.testing.assert_allclose(
                np.asarray(art.axes[index].collections[0].get_edgecolor(), dtype=float)[
                    0
                ],
                mpl.colors.to_rgba(CONVERGENCE_COLORS[mesh]),
            )
            assert mesh in art.axes[index].get_title()
        for index, mesh in enumerate(CONVERGENCE_DISPLAY_MESHES):
            np.testing.assert_allclose(
                np.asarray(art.axes[6].collections[index].get_edgecolor(), dtype=float)[
                    0
                ],
                mpl.colors.to_rgba(CONVERGENCE_COLORS[mesh]),
            )
        assert len(set(CONVERGENCE_COLORS.values())) == 6
        assert all("50000" not in axis.get_title() for axis in art.axes)
        assert art.axes[17].get_yscale() == "linear"
        assert art.axes[22].get_yscale() == "symlog"
        assert art.axes[-1].get_yscale() == "linear"
        assert any("no extrapolation" in text.get_text() for text in art.texts)

    monkeypatch.setattr(
        "research.spt.publication_timecourses.save_publication", inspect
    )
    render_convergence(data, sources, tmp_path, RenderOptions(formats=("svg",), dpi=30))


def test_cached_publication_skips_arrays_and_repairs_missing_exports(
    publication_data: PublicationData,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Unchanged figure generation uses metadata only; missing PNGs are regenerated."""
    data = publication_data
    source = data.source("simulations_lobulus", "sim014")
    monkeypatch.setattr(
        "research.spt.publication.simulation_cohort", lambda data: [source]
    )
    first = create_publication(data.root, figures=("4",), formats=("png",), dpi=30)
    pointer = data.root / "publication/figure-04/manifest.json"
    before = pointer.read_bytes()
    exported = data.root / "publication/exports/figure-04.png"
    assert exported.samefile(first["4"] / "figure-04.png")
    exported.unlink()
    with monkeypatch.context() as cached:
        cached.setattr(
            "numpy.load",
            lambda *args, **kwargs: pytest.fail("Cache hit opened mesh arrays"),
        )
        cached.setattr(
            "xarray.open_dataset",
            lambda *args, **kwargs: pytest.fail("Cache hit opened numerical arrays"),
        )
        assert (
            create_publication(data.root, figures=("4",), formats=("png",), dpi=30)
            == first
        )
    assert pointer.read_bytes() == before
    assert exported.samefile(first["4"] / "figure-04.png")
    (first["4"] / "figure-04.png").unlink()
    again = create_publication(data.root, figures=("4",), formats=("png",), dpi=30)
    assert (again["4"] / "figure-04.png").is_file()
    provenance = json.loads((again["4"] / "provenance.json").read_text())
    assert provenance["snapshot_times_s"] == list(SIMULATION_TIMES)
    assert provenance["conditions"][0]["sim_id"] == "sim014"
    assert provenance["raw_hdf_reads"] == 0


def test_compact_exports_preserve_full_size_and_resolution(
    publication_data: PublicationData, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Compact exports retain 300 DPI, original scalar styles and independent caches."""
    data = publication_data
    source = data.source("simulations_lobulus", "sim014")
    monkeypatch.setattr(
        "research.spt.publication.simulation_cohort", lambda data: [source]
    )
    original = create_publication(data.root, figures=("4",), formats=("png",), dpi=30)
    pointer = data.root / "publication/figure-04/manifest.json"
    before = pointer.read_bytes()
    compact = create_publication(
        data.root, figures=("4",), formats=("png",), dpi=300, width_mm=180
    )
    exported = data.root / "publication/exports/figure-04-compact.png"
    assert exported.samefile(compact["4"] / exported.name)
    assert pointer.read_bytes() == before
    assert original["4"] != compact["4"]
    with Image.open(exported) as image:
        assert image.info["dpi"] == pytest.approx((300, 300), abs=0.01)
        assert image.width == pytest.approx(180 / 25.4 * 300, rel=0.05)
    provenance = json.loads((compact["4"] / "provenance.json").read_text())
    assert provenance["width_mm"] == 180 and provenance["variant"] == "compact"
    assert provenance["snapshot_times_s"] == list(SIMULATION_TIMES)
    with monkeypatch.context() as cached:
        cached.setattr("numpy.load", _no_hdf)
        cached.setattr("xarray.open_dataset", _no_hdf)
        assert (
            create_publication(
                data.root, figures=("4",), formats=("png",), dpi=300, width_mm=180
            )
            == compact
        )


def test_mixed_readouts_sample_native_cells_at_one_common_time(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Unequal sector meshes retain cell order, interpolation and discrete states."""
    raw, root = tmp_path / "raw", tmp_path / "local"
    for index, mesh in enumerate(CONVERGENCE_MESHES):
        _input(raw, f"convergence_{mesh}", index + 1, (0, 1, 2))
    process_geometries(raw, root, study="convergence", workers=1)
    process_studies(raw, root, study="convergence", geometry="lobulus", workers=1)
    sources = plot_sources(raw, root, study="convergence", geometry="lobulus")
    for source in sources:
        prepare_analysis(
            source, root, root, AnalysisOptions(time_points=3), (0, 1, 2), (0, 2)
        )
        prepare_samples(source, root, AnalysisOptions().fields, (0, 2))
    mixed = plot_sources(raw, root, study="convergence", geometry="combined")[0]
    prepare_samples(mixed, root, AnalysisOptions().fields, (0, 2))
    data = PublicationData(root)
    monkeypatch.setattr("h5py.File", _no_hdf)
    monkeypatch.setattr("meshio.xdmf.TimeSeriesReader", _no_hdf)
    result = data.mixed_snapshot(0.5)
    np.testing.assert_array_equal(result.times, [0.5])
    for sector in range(6):
        selected = result.sectors == sector
        assert np.count_nonzero(selected) == sector + 1
        np.testing.assert_allclose(result.values[0, 0, selected], sector + 2)
        np.testing.assert_allclose(
            result.values[0, 2, selected], np.linspace(0.2, 0.8, sector + 1)
        )
        np.testing.assert_array_equal(result.values[0, -1, selected], 0)
    displayed = convergence_view(result)
    assert displayed.cell_ids is result.cell_ids
    assert displayed.values is result.values
    assert displayed.times is result.times
    # Resolutions and their scalar values move together, without changing caches.
    np.testing.assert_array_equal(data.mixed_snapshot(0.5).faces, result.faces)
    for sector, mesh in enumerate(CONVERGENCE_DISPLAY_MESHES):
        original = CONVERGENCE_MESHES.index(mesh)
        selected = displayed.sectors == sector
        assert np.count_nonzero(selected) == original + 1
        np.testing.assert_allclose(displayed.values[0, 0, selected], original + 2)
        centers = displayed.faces[displayed.sectors[displayed.cell_ids] == sector].mean(
            axis=(0, 1)
        )
        angle = np.rad2deg(np.arctan2(centers[1], centers[0])) % 360
        assert angle == pytest.approx((120 - 60 * sector) % 360)
    np.testing.assert_array_equal(data.mixed_snapshot(2).values[0, -1], 1)
    with pytest.raises(ValueError, match="outside native history"):
        data.mixed_snapshot(3)
    geometry_pointer = mixed.run.xdmf.parent.parent / "manifest.json"
    composition = json.loads(geometry_pointer.read_text())
    composition["configuration"]["sources"].reverse()
    geometry_pointer.write_text(json.dumps(composition))
    # The geometry sample validation catches the changed local composition first.
    with pytest.raises(ValueError, match=r"stale|source order"):
        data.mixed_snapshot(0.5)


def test_publication_figures_appear_in_fast_report(tmp_path: Path) -> None:
    """Expose publication figures in their own report section with small previews."""
    (tmp_path / "catalog.tsv").write_text(
        "run_id\tpattern_name\tboundary_flow_key\n"
        "sim001\tconstant\t0\nsim014\tlinear increase\t5\n"
    )
    for number in (7, 6, 5, 4):
        stem = f"figure-{number:02}"
        _publish(
            tmp_path,
            f"publication/{stem}",
            (f"{stem}.png", f"{stem}.pdf", "provenance.json"),
            {
                "study": "simulations_lobulus",
                "title": f"Figure {number}",
                "runs": ["sim014"] if number == 4 else ["sim001", "sim014"],
            },
        )
    _publish(
        tmp_path,
        "publication/figure-04-compact",
        ("figure-04-compact.png",),
        {
            "study": "simulations_lobulus",
            "title": "Figure 4 (compact)",
            "runs": ["sim014"],
        },
    )
    report = ReportParser(build_report(tmp_path).read_text())
    assert [card["data-title"] for card in report.cards] == [
        "Figure 4",
        "Figure 4 (compact)",
        "Figure 5",
        "Figure 6",
        "Figure 7",
    ]
    assert all(card["data-section"] == "publication" for card in report.cards)
    assert report.images[0]["src"] and str(report.images[0]["src"]).endswith(".webp")


@pytest.mark.parametrize(
    "figures,dose", [((), 6), (("4", "4"), 6), (("8",), 6), (("4",), 0), (("4",), 9)]
)
def test_invalid_selection_fails_before_output(
    tmp_path: Path, figures: tuple[str, ...], dose: int
) -> None:
    """Reject invalid figure identifiers and doses before any output is created."""
    with pytest.raises(ValueError):
        create_publication(tmp_path, figures=figures, dose=dose)
    assert not list(tmp_path.iterdir())


def test_content_crop_keeps_white_tissue_and_annotations(
    publication_data: PublicationData, tmp_path: Path
) -> None:
    """Crop empty axis slots while preserving white tissue, labels and vector text."""
    from research.spt.fields import LAYERS
    from research.spt.plot_render import figure

    source = publication_data.source("simulations_lobulus", "sim014")
    samples = publication_data.snapshots(source, (0.0,))
    art = figure((8, 6))
    art.set_layout_engine("none")
    try:
        axis = art.add_axes((0.3, 0.3, 0.4, 0.4))
        raster_panel(
            axis,
            samples,
            replace(LAYERS[-1], show_edges=False),
            samples.values[0, -1],
            (0, 1),
            side=2,
            dpi=100,
        )
        pixels = axis.images[0].get_array()
        assert pixels is not None and np.all(pixels[..., :3] == 255)
        axis.set_title("White living tissue")
        axis.annotate(
            "Retain this label",
            xy=(0, 0.5),
            xytext=(-0.4, 0.5),
            xycoords="axes fraction",
            annotation_clip=False,
            arrowprops={"arrowstyle": "->"},
        )
        with mpl.rc_context({"svg.fonttype": "none"}):
            save_publication(
                art,
                tmp_path / "cropped",
                RenderOptions(formats=("png", "svg"), dpi=100),
            )
        with Image.open(tmp_path / "cropped.png") as exported:
            assert exported.width < 600 and exported.height < 400
            assert exported.height > 200
            assert exported.info["dpi"] == pytest.approx((100, 100), abs=0.01)
        svg = (tmp_path / "cropped.svg").read_text()
        assert "Retain this label" in svg and "White living tissue" in svg
        assert "<text" in svg
    finally:
        art.clear()


def test_retire_only_obsolete_compact_exports(tmp_path: Path) -> None:
    """Remove obsolete compact cards and exports without touching full-size figures."""
    for stem, number, variant in (
        ("figure-06-compact", "6", "compact"),
        ("figure-04-compact", "4", "compact"),
        ("figure-06", "6", "full"),
    ):
        generation = _publish(
            tmp_path,
            f"publication/{stem}",
            (f"{stem}.png",),
            {"figure": number, "variant": variant, "study": "simulations_lobulus"},
        )
        directory = tmp_path / "publication/exports"
        directory.mkdir(exist_ok=True)
        (directory / f"{stem}.png").hardlink_to(generation / f"{stem}.png")
    remove_obsolete_compact(tmp_path)
    assert not (tmp_path / "publication/figure-06-compact/manifest.json").exists()
    assert not (tmp_path / "publication/exports/figure-06-compact.png").exists()
    assert (tmp_path / "publication/exports/figure-06.png").is_file()
    assert (tmp_path / "publication/exports/figure-04-compact.png").is_file()
    remove_obsolete_compact(tmp_path)


def test_crop_excludes_padding_of_frameless_legends(tmp_path: Path) -> None:
    """The crop reaches the legend title rather than its invisible outer frame."""
    from research.spt.plot_render import figure

    art = figure((6, 4))
    try:
        axis = art.subplots()
        (line,) = axis.plot([0, 1], [0, 1], label="06250")
        art.legend(
            handles=[line],
            title="Resolution",
            loc="outside upper center",
            frameon=False,
        )
        save_publication(
            art, tmp_path / "legend", RenderOptions(formats=("png",), dpi=100)
        )
        with Image.open(tmp_path / "legend.png") as image:
            pixels = np.asarray(image.convert("RGB"))
            occupied = np.flatnonzero((pixels < 255).any(axis=(1, 2)))
            assert occupied[0] <= 5
    finally:
        art.clear()


@pytest.mark.parametrize("figure", ["6", "7", "S1"])
def test_compact_selection_rejects_other_figures(tmp_path: Path, figure: str) -> None:
    """Only the large geometry composites need a compact publication variant."""
    with pytest.raises(ValueError, match="only for Figures 4 and 5"):
        create_publication(tmp_path, figures=(figure,), width_mm=180)
    assert not list(tmp_path.iterdir())


@pytest.mark.parametrize(
    "width,time",
    [(0, 5400), (-1, 5400), (float("nan"), 5400), (180, -1), (180, float("inf"))],
)
def test_invalid_dimensions_or_common_time_fail_before_output(
    tmp_path: Path, width: float, time: float
) -> None:
    """Reject invalid physical dimensions and times without touching existing caches."""
    with pytest.raises(ValueError):
        create_publication(tmp_path, width_mm=width, comparison_time=time)
    assert not list(tmp_path.iterdir())
