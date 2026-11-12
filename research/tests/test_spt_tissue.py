"""Verify complete square coverage, native cell identity and flat tissue exports."""

import json
from collections.abc import Callable
from dataclasses import replace
from pathlib import Path
from unittest.mock import patch

import matplotlib as mpl
import numpy as np
import pytest
from PIL import Image

from research.spt.fields import LAYERS
from research.spt.manuscript_data import ManuscriptData
from research.spt.report import build_report
from research.spt.tissue_geometry import clip_square, crop_tissue, tile_lobules
from research.spt.tissue_panels import (
    TissueOptions,
    create_tissue_panels,
    main,
    panel_stem,
    patch_image,
)
from research.tests.test_spt_manuscript import manuscript_data as manuscript_data
from research.tests.test_spt_report import ReportParser


def _hexagon(hole: bool = False) -> tuple[np.ndarray, np.ndarray]:
    angles = np.pi / 6 + np.arange(6) * np.pi / 3
    vertices = np.column_stack((np.cos(angles), np.sin(angles)))
    faces = []
    for index in range(6):
        first, second = vertices[index], vertices[(index + 1) % 6]
        if hole:
            faces.extend(
                [[first, second, first * 0.1], [second, second * 0.1, first * 0.1]]
            )
        else:
            faces.append([[0, 0], first, second])
    return np.array(faces), np.arange(len(faces))


@pytest.mark.parametrize("width, lobules", [(2.0, 7), (2.5, 17)])
def test_full_lattice_covers_square_and_preserves_partial_cells(
    width: float, lobules: int
) -> None:
    """Cover the whole square and its corners, with exact cut edges and native IDs."""
    faces, ids = _hexagon()
    before = faces.copy()
    patch = tile_lobules(faces, ids, width)
    half = width / 2
    assert len(patch.centers) == lobules
    assert patch.width == width
    assert patch.source_width == pytest.approx(np.sqrt(3))
    assert patch.nominal_area == pytest.approx(width**2)
    assert patch.tissue_area == pytest.approx(width**2)
    np.testing.assert_array_equal(patch.faces.min(axis=(0, 1)), [-half, -half])
    np.testing.assert_array_equal(patch.faces.max(axis=(0, 1)), [half, half])
    assert set(patch.cell_ids) == set(ids)
    for lobule in np.unique(patch.lobule_ids):
        selected = patch.lobule_ids == lobule
        assert np.all(patch.faces[selected] <= half)
        assert np.all(patch.faces[selected] >= -half)
    native = faces[patch.cell_ids] / patch.source_width
    local = patch.faces - patch.centers[patch.lobule_ids, None]
    edges = native[:, 1:] - native[:, :1]
    weights = np.linalg.solve(
        np.swapaxes(edges, 1, 2), np.swapaxes(local - native[:, :1], 1, 2)
    )
    assert np.all(weights >= -1e-12)
    assert np.all(weights.sum(axis=1) <= 1 + 1e-12)
    np.testing.assert_array_equal(faces, before)
    clipped = clip_square(np.array([[-2, -2], [2, -2], [2, 2], [-2, 2]]), half)
    assert len(clipped) == 4
    np.testing.assert_array_equal(np.abs(clipped), np.full((4, 2), half))


def test_vessel_lumens_are_retained_and_never_assigned_tissue_values() -> None:
    """Periodic tissue closes outer margins while the native central lumen stays open."""
    faces, ids = _hexagon(hole=True)
    patch = tile_lobules(faces, ids)
    assert patch.tissue_area < patch.nominal_area
    options = TissueOptions(dpi=40, vessel_color="lightgray")
    pixels = patch_image(patch, np.ones(len(ids)), LAYERS[-1], (0, 1), options)
    np.testing.assert_array_equal(pixels[60, 60], [211, 211, 211, 255])
    assert np.any(np.all(pixels[..., :3] == 0, axis=-1))
    assert np.all(pixels[..., 3] == 255)


@pytest.mark.parametrize("center", [(0.0, 0.0), (0.371, -0.482), (-1.0, 1.0)])
def test_offset_crop_preserves_coverage_and_native_cell_identity(
    center: tuple[float, float],
) -> None:
    """Crop a full square from larger tissue without moving or averaging cells."""
    faces, ids = _hexagon()
    tissue = tile_lobules(faces, ids, 4)
    before = tissue.faces.copy()
    cropped = crop_tissue(tissue, 2, center)
    assert cropped.width == 2
    assert cropped.nominal_area == pytest.approx(4)
    assert cropped.tissue_area == pytest.approx(4)
    np.testing.assert_array_equal(cropped.faces.min(axis=(0, 1)), [-1, -1])
    np.testing.assert_array_equal(cropped.faces.max(axis=(0, 1)), [1, 1])
    native = faces[cropped.cell_ids] / cropped.source_width
    local = cropped.faces - cropped.centers[cropped.lobule_ids, None]
    edges = native[:, 1:] - native[:, :1]
    weights = np.linalg.solve(
        np.swapaxes(edges, 1, 2), np.swapaxes(local - native[:, :1], 1, 2)
    )
    assert np.all(weights >= -1e-12)
    assert np.all(weights.sum(axis=1) <= 1 + 1e-12)
    np.testing.assert_array_equal(tissue.faces, before)
    pixels = patch_image(
        cropped, np.ones(len(ids)), LAYERS[-1], (0, 1), TissueOptions(dpi=40)
    )
    assert np.all(pixels[..., :3] == 0)
    assert np.all(pixels[..., 3] == 255)


def test_offset_crop_keeps_native_vessel_lumens() -> None:
    """Shift the central vessel with the tissue, preserving its empty lumen."""
    faces, ids = _hexagon(hole=True)
    tissue = tile_lobules(faces, ids, 4)
    cropped = crop_tissue(tissue, 2, (0.3, -0.4))
    assert cropped.tissue_area < cropped.nominal_area
    pixels = patch_image(
        cropped,
        np.ones(len(ids)),
        LAYERS[-1],
        (0, 1),
        TissueOptions(dpi=40, vessel_color="lightgray"),
    )
    np.testing.assert_array_equal(pixels[36, 42], [211, 211, 211, 255])


@pytest.mark.parametrize(
    "center", [(float("nan"), 0.0), (0.0, float("inf")), (1.01, 0.0)]
)
def test_invalid_or_outside_crop_fails(center: tuple[float, float]) -> None:
    """Reject crops that cannot preserve a complete, finite tissue square."""
    faces, ids = _hexagon()
    tissue = tile_lobules(faces, ids, 4)
    with pytest.raises(ValueError, match=r"center|inside"):
        crop_tissue(tissue, 2, center)
    with pytest.raises(ValueError, match="inside"):
        crop_tissue(tissue, 5, (0, 0))


@pytest.mark.parametrize("width", [2.0, 2.5])
def test_flat_square_pixels_cover_every_corner_and_repeat_native_colors(
    width: float,
) -> None:
    """Top-down normalized coordinates map directly to pixels, without camera skew."""
    faces, ids = _hexagon()
    patch = tile_lobules(faces, ids, width)
    half = width / 2
    options = TissueOptions(dpi=80)
    black = patch_image(patch, np.ones(len(ids)), LAYERS[-1], (0, 1), options)
    assert black.shape == (240, 240, 4)
    assert np.all(black[..., :3] == 0)
    assert np.all(black[..., 3] == 255)
    values = np.linspace(0, 1, len(ids))
    pixels = patch_image(patch, values, LAYERS[2], (0, 1), options)
    for x in (-1, 0, 1):
        for face in range(6):
            point = faces[face].mean(axis=0) / patch.source_width + [x, 0]
            if np.max(np.abs(point)) >= half:
                continue
            column = int((point[0] + half) / width * 240)
            row = int((half - point[1]) / width * 240)
            expected = np.array(mpl.colormaps[LAYERS[2].colormap](values[face])) * 255
            np.testing.assert_allclose(pixels[row, column], expected, atol=1)


@pytest.mark.parametrize("width", [0.0, 1.99, float("nan"), float("inf")])
def test_invalid_square_width_fails(width: float) -> None:
    """Reject undersized or nonfinite domains before building the lattice."""
    faces, ids = _hexagon()
    with pytest.raises(ValueError, match="at least 2 lobules"):
        tile_lobules(faces, ids, width)


def test_nonhexagonal_or_corrupt_geometry_is_rejected() -> None:
    """Reject malformed source cells and meshes that cannot tile without gaps."""
    faces, ids = _hexagon()
    with pytest.raises(ValueError, match=r"hexagonal|regular"):
        tile_lobules(faces * [1.1, 1], ids)
    with pytest.raises(ValueError, match="cell ID"):
        tile_lobules(faces, ids - 1)
    faces[0, 0, 0] = np.nan
    with pytest.raises(ValueError, match="finite"):
        tile_lobules(faces, ids)


def test_fractional_time_filenames_are_unique_and_keep_their_extension() -> None:
    """Keep distinct fractional times in separate correctly named image files."""
    first = panel_stem("tissue", LAYERS[2], 0.5)
    second = panel_stem("tissue", LAYERS[2], 0.500000001)
    assert first != second
    assert Path(first).suffix == Path(second).suffix == ""
    assert Path(first).with_suffix(".png").stem == first


@pytest.mark.parametrize("vessel_color", ("white", "lightgray"))
def test_local_tissue_workflow_exports_and_cached_reads(
    manuscript_data: ManuscriptData,
    monkeypatch: pytest.MonkeyPatch,
    reject_history_hdf: Callable[[str], None],
    vessel_color: str,
) -> None:
    """Export real synthetic native cells, exact square files and a reusable report."""
    data = manuscript_data
    source = data.source("simulations_lobulus", "sim014")
    monkeypatch.setattr(
        "research.spt.tissue_panels.simulation_cohort", lambda data: [source]
    )
    options = TissueOptions(times=(0, 86400), dpi=40, vessel_color=vessel_color)

    def forbidden(*args: object, **kwargs: object) -> None:
        pytest.fail("Tissue panels must not reopen cached arrays")

    reject_history_hdf("Tissue panels must not reopen mesh or HDF histories")
    generation = create_tissue_panels(data.root, options)
    pngs = sorted(generation.glob("*.png"))
    assert len(pngs) == 7
    for path in pngs:
        with Image.open(path) as image:
            if path.stem == "linear-increase-tissue-d6":
                continue
            assert image.size == (120, 120)
            assert image.info["dpi"] == pytest.approx((40, 40), abs=0.01)
    provenance = json.loads((generation / "provenance.json").read_text())
    assert provenance["geometry_view"]["tilt_from_top_degrees"] == 0
    assert provenance["options"]["width"] == 2.0
    assert provenance["options"]["seed"] == 0
    assert provenance["layout"]["tissue_width_in_lobule_widths"] == 4
    layouts = provenance["layout"]["panels"]
    assert set(layouts) == {
        path.stem for path in pngs if path.stem != "linear-increase-tissue-d6"
    }
    centers = [
        tuple(layout["crop_center_in_lobule_widths"]) for layout in layouts.values()
    ]
    assert len(set(centers)) == 6
    assert np.all(np.abs(centers) <= 1)
    for layout in layouts.values():
        assert layout["square_area_in_lobule_widths_squared"] == pytest.approx(4)
        assert layout["clipped_triangles"] > 0
    assert provenance["options"]["times"] == [0, 86400]
    assert provenance["scalar_ranges"]["rr_necrosis"] == [0, 1]
    assert provenance["raw_hdf_reads"] == 0
    svg = (generation / "linear-increase-tissue-d6.svg").read_text()
    assert ("vessel lumen" in svg) == (vessel_color == "lightgray")
    document = build_report(data.root).read_text()
    assert "2 x 2 lobule widths" in document
    parsed = ReportParser(document)
    assert len(parsed.cards) == 7
    assert "CYP2E1 at 24 h" in document and "Necrosis at 24 h" in document
    assert document.index("CYP2E1 at 0 h") < document.index("CYP2E1 at 24 h")
    assert document.index("CYP2E1 at 24 h") < document.index(
        "Toxic compound T [mM] at 0 h"
    )
    titles = [provenance["artifact_titles"][p.stem] for p in pngs]
    assert len(set(titles)) == 7
    target = (
        data.root / "manuscript/exports/linear-increase-tissue-d6-cyp2e1-t86400s.png"
    )
    target.unlink()
    with patch("numpy.load", forbidden), patch("xarray.open_dataset", forbidden):
        assert create_tissue_panels(data.root, options) == generation
    assert target.is_file()


def test_random_crop_seed_reproduces_positions_and_changes_exports(
    manuscript_data: ManuscriptData, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Reproduce random crops with the same seed and sample new crops with another."""
    data = manuscript_data
    source = data.source("simulations_lobulus", "sim014")
    monkeypatch.setattr(
        "research.spt.tissue_panels.simulation_cohort", lambda data: [source]
    )
    options = TissueOptions(times=(0, 86400), dpi=40, formats=("png",))
    first = create_tissue_panels(data.root, options)
    first_layout = json.loads((first / "provenance.json").read_text())["layout"]
    first_images = {
        name: (first / f"{name}.png").read_bytes() for name in first_layout["panels"]
    }
    second = create_tissue_panels(data.root, options, force=True)
    assert second != first
    second_layout = json.loads((second / "provenance.json").read_text())["layout"]
    assert second_layout == first_layout
    for name, pixels in first_images.items():
        assert pixels == (second / f"{name}.png").read_bytes()
    changed = create_tissue_panels(data.root, replace(options, seed=1))
    changed_layout = json.loads((changed / "provenance.json").read_text())["layout"]
    assert changed_layout != first_layout
    assert any(
        pixels != (changed / f"{name}.png").read_bytes()
        for name, pixels in first_images.items()
    )


def test_invalid_options_fail_before_writing(tmp_path: Path) -> None:
    """Validate selections before creating any generation or output directory."""
    for options in (
        TissueOptions(dose=0),
        TissueOptions(width=1.99),
        TissueOptions(seed=-1),
        TissueOptions(times=(1, 0)),
        TissueOptions(times=(float("nan"),)),
        TissueOptions(vessel_color="red"),
    ):
        with pytest.raises(ValueError):
            create_tissue_panels(tmp_path, options)
    assert not list(tmp_path.iterdir())
    assert main(["--data-root", str(tmp_path), "--width-lobules", "1.99"]) == 1
    assert main(["--data-root", str(tmp_path), "--seed", "-1"]) == 1
    assert not list(tmp_path.iterdir())
