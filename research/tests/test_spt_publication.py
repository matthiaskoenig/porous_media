"""Check publication annotations, true sampling times, local reads and cache reuse."""

import json
import xml.etree.ElementTree as ET
from dataclasses import replace
from pathlib import Path

import matplotlib as mpl
import numpy as np
import pytest
import xarray as xr
from matplotlib.figure import Figure

from research.spt.analysis_data import AnalysisOptions, prepare_analysis
from research.spt.conditions import simulation_conditions
from research.spt.data_processing import process_studies
from research.spt.geometry import process_geometries
from research.spt.plot_data import PlotSource, plot_sources, prepare_samples
from research.spt.plot_render import RenderOptions, Samples
from research.spt.plot_style import CONVERGENCE_ORDER
from research.spt.publication import create_publication, simulation_cohort
from research.spt.publication_data import SIMULATION_TIMES, PublicationData
from research.spt.publication_geometry import (
    raster_panel,
    render_example,
    render_necrosis,
    render_patterns,
)
from research.spt.publication_timecourses import render_convergence
from research.spt.report import build_report
from research.tests.test_spt_analysis import _input
from research.tests.test_spt_report import ReportParser, _publish


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
        assert pixels.shape == (150, 150, 4)
        assert pixels.nbytes == 150 * 150 * 4
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
        label in patterns
        for label in [
            "A  constant",
            "B  linear increase",
            "C  sharp perivenous",
            "D  linear decrease",
            "E  sharp periportal",
            "F  random",
        ]
    )
    necrosis = (tmp_path / "figure-07.svg").read_text()
    assert all(f"Dose {dose}" in necrosis for dose in range(1, 9))
    assert "Substrate / drug inflow" in necrosis


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
        for mesh in CONVERGENCE_ORDER
    ]
    monkeypatch.setattr(data, "mesh", lambda study, run_id: samples)

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

    def inspect(art: Figure, destination: Path, options: RenderOptions) -> None:
        assert destination.name == "figure-s1-convergence"
        assert len(art.axes) == 18
        assert len(art.axes[0].images) == 1
        for row in (8, 13):
            curves = {
                line.get_label(): line
                for line in art.axes[row].lines
                if line.get_label() in CONVERGENCE_ORDER
            }
            assert len(curves) == 7
            assert np.asarray(curves["06250"].get_xdata())[-1] == 36590 / 3600
            assert np.asarray(curves["09375"].get_xdata())[-1] == 77687.6 / 3600
            assert np.asarray(curves["50000"].get_xdata())[-1] == 24
        assert art.axes[8].get_yscale() == "linear"
        assert art.axes[13].get_yscale() == "symlog"
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
    report = ReportParser(build_report(tmp_path).read_text())
    assert [card["data-title"] for card in report.cards] == [
        f"Figure {n}" for n in (4, 5, 6, 7)
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
