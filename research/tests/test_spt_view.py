"""Check native thickness, scientific colors, exterior walls and cache adoption."""

from dataclasses import replace
from pathlib import Path
from unittest.mock import patch

import meshio
import numpy as np
import pytest

from research.spt.artifacts import generation_path
from research.spt.fields import LAYERS
from research.spt.plot_data import load_heights, plot_sources, prepare_samples
from research.spt.plot_mesh import boundary_edges, exterior_edges
from research.spt.plot_render import Samples, figure, panel
from research.spt.plot_view import (
    SurfaceProjection,
    camera_matrix,
    face_heights,
    project,
)
from research.tests.test_spt_geometry import _source


@pytest.mark.parametrize("cell_type", ("hexahedron", "wedge"))
def test_tilt_keeps_native_thickness_and_cell_colors(cell_type: str) -> None:
    """Caps and walls keep their owning cells, canonical colors and true thickness."""
    bottom = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]], dtype=float)
    points = np.vstack((bottom, bottom + np.array([0, 0, 0.2])))
    faces = bottom[[[0, 1, 3], [1, 2, 3]], :2]
    connectivity = np.array([[0, 1, 3, 4, 5, 7], [1, 2, 3, 5, 6, 7]])
    cells = [
        meshio.CellBlock(
            cell_type, [np.arange(8)] if cell_type == "hexahedron" else connectivity
        )
    ]
    cell_ids = np.array([0, 0] if cell_type == "hexahedron" else [0, 1])
    values = np.array([0.4] if cell_type == "hexahedron" else [0.2, 0.8])
    heights = face_heights(points, cells)
    np.testing.assert_array_equal(heights, [[0, 0.2], [0, 0.2]])
    np.testing.assert_allclose(
        camera_matrix() @ camera_matrix().T, np.eye(3), atol=1e-15
    )
    data = Samples(
        faces, cell_ids, np.empty(0), np.zeros(len(values)), np.empty(0), heights
    )
    caps = data.surface.caps
    np.testing.assert_allclose(caps, project(faces, 0.2))
    np.testing.assert_allclose(
        caps - project(faces, 0),
        np.broadcast_to(0.2 * camera_matrix()[:2, 2], caps.shape),
    )
    assert np.linalg.norm(caps[0, 1] - caps[0, 0]) < 1
    assert len(data.surface.walls[0]) == 2
    assert len(data.display_ids) == 4
    assert set(data.display_ids[:2]) == set(cell_ids)
    np.testing.assert_array_equal(data.display_ids[2:], cell_ids)
    art = figure((3, 3))
    try:
        axis = art.subplots()
        layer = replace(LAYERS[2], show_edges=False)
        collection = panel(axis, data, layer, values, (0, 1))
        art.canvas.draw()
        np.testing.assert_array_equal(collection.get_array(), values[data.display_ids])
        np.testing.assert_allclose(
            np.asarray(collection.get_facecolor(), dtype=float),
            collection.cmap(values[data.display_ids]),
        )
        np.testing.assert_array_equal(data.faces, faces)
        np.testing.assert_array_equal(data.heights, heights)
        painted = np.concatenate(data.surface.polygons)
        low, high = data.surface.bounds
        assert np.all(painted >= low) and np.all(painted <= high)
        assert axis.get_xlim()[0] <= low[0] and axis.get_xlim()[1] >= high[0]
        assert axis.get_ylim()[0] <= low[1] and axis.get_ylim()[1] >= high[1]
    finally:
        art.clear()


def test_holes_keep_real_walls_without_internal_cell_walls() -> None:
    """A triangulated ring exposes its hole and exterior, never internal diagonals."""
    faces = []
    for x in range(3):
        for y in range(3):
            if x == y == 1:
                continue
            quad = np.array(
                [[x, y], [x + 1, y], [x + 1, y + 1], [x, y + 1]], dtype=float
            )
            faces.extend(quad[[[0, 1, 2], [0, 2, 3]]])
    native = np.array(faces)
    edges, owners = boundary_edges(native)
    assert len(edges) == 16
    for edge, owner in zip(edges, owners, strict=True):
        assert all(
            any(np.array_equal(vertex, point) for point in native[owner])
            for vertex in edge
        )
    view = SurfaceProjection(native, np.tile([0, 0.2], (len(native), 1)))
    walls, wall_owners = view.walls
    assert len(walls) == 8
    assert set(wall_owners).issubset(owners)
    midpoints = native[wall_owners].mean(axis=1)
    hole_walls = (
        (midpoints[:, 0] < 1) & (midpoints[:, 1] > 1) & (midpoints[:, 1] < 2)
    ) | ((midpoints[:, 0] > 1) & (midpoints[:, 0] < 2) & (midpoints[:, 1] > 2))
    assert np.count_nonzero(hole_walls) == 2
    assert len(view.face_indices) == len(native) + len(walls)


def test_mixed_resolution_seams_never_become_internal_walls() -> None:
    """Unequal edge divisions join as solid tissue; exposed partial edges remain."""

    def triangles(x: float, bottom: float, top: float) -> np.ndarray:
        quad = np.array([[x, bottom], [x + 0.5, bottom], [x + 0.5, top], [x, top]])
        return quad[[[0, 1, 2], [0, 2, 3]]]

    left = triangles(0, 0, 1)
    lower_right = triangles(0.5, 0, 0.5)
    upper_right = triangles(0.5, 0.5, 1)
    faces = np.concatenate((left, lower_right, upper_right))
    assert len(boundary_edges(faces)[0]) == 10
    edges, owners = exterior_edges(faces)
    assert len(edges) == 7 and len(owners) == 7
    assert not np.any(np.all(edges[:, :, 0] == 0.5, axis=1))
    view = SurfaceProjection(faces, np.tile([0, 0.2], (len(faces), 1)))
    assert len(view.walls[0]) == 4

    partial = np.concatenate((left, lower_right))
    edges, _ = exterior_edges(partial)
    exposed = edges[np.all(edges[:, :, 0] == 0.5, axis=1)]
    assert len(exposed) == 1
    np.testing.assert_array_equal(np.sort(exposed[0, :, 1]), [0.5, 1])


def test_old_samples_gain_depth_without_resampling_scalars(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Adopt an old planar NPZ using geometry only, then resume without mesh reads."""
    raw, output = tmp_path / "raw", tmp_path / "output"
    descriptor = _source(raw, "sim001", 3, (0, 2, 4))
    source = plot_sources(raw, output, study="simulations", geometry="sixth")[0]
    manifest = prepare_samples(source, output, ("rr_(T)",), (0, 2, 4))
    assert manifest is not None
    original = (
        generation_path(source.directory(output, "plot-data"), manifest) / "samples.npz"
    )
    legacy = tmp_path / "planar.npz"
    with np.load(original, allow_pickle=False) as data:
        np.savez(legacy, **{key: data[key] for key in data.files if key != "heights"})
    before = {
        path: (path.stat().st_size, path.stat().st_mtime_ns)
        for path in (descriptor, original, legacy)
    }
    directory = source.directory(output, "plot-depth")
    with patch.object(
        meshio.xdmf.TimeSeriesReader,
        "read_data",
        side_effect=AssertionError("Scalar histories must stay closed"),
    ):
        heights = load_heights(directory, manifest, raw, sample_path=legacy)
    np.testing.assert_array_equal(heights, np.tile([0, 0.2], (3, 1)))
    monkeypatch.setattr(
        "meshio.xdmf.TimeSeriesReader",
        lambda *args: pytest.fail("Cached depth reopened a mesh"),
    )
    np.testing.assert_array_equal(
        load_heights(directory, manifest, raw, sample_path=legacy), heights
    )
    np.testing.assert_array_equal(
        load_heights(directory, manifest, raw, sample_path=original), heights
    )
    assert before == {
        path: (path.stat().st_size, path.stat().st_mtime_ns) for path in before
    }
