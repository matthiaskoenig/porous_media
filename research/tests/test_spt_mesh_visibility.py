"""Check that dense mesh grids do not obscure convergence scalar fields."""

from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest
from matplotlib.collections import LineCollection
from matplotlib.figure import Figure

from research.spt.fields import LAYERS
from research.spt.plot_mesh import boundary_segments
from research.spt.plot_render import (
    RenderOptions,
    Samples,
    figure,
    panel,
    render_mesh,
    render_timecourse,
)
from research.tests.test_spt_plot_style import _source


def _square() -> Samples:
    return Samples(
        np.array([[[0, 0], [1, 0], [0, 1]], [[1, 0], [1, 1], [0, 1]]], dtype=float),
        np.array([0, 1]),
        np.array([[[0.2, 0.8]]]),
        np.array([0, 1]),
        np.array([0]),
    )


@pytest.mark.parametrize("offset", [1e-13, 1e-6])
def test_boundary_edges_exclude_interior_diagonals_and_join_seams(
    offset: float,
) -> None:
    """Coordinate matching removes artificial seams despite reflection roundoff."""
    samples = _square()
    samples.faces[1] += offset
    assert len(samples.boundary) == 4
    assert len(samples.sector_boundaries) == 6
    assert samples.boundary is samples.boundary
    assert not any(np.allclose(edge, [[1, 0], [0, 1]]) for edge in samples.boundary)


def test_boundary_extraction_preserves_holes() -> None:
    """Keep the inner hole as well as the external perimeter of a ring mesh."""
    faces = []
    for x in range(3):
        for y in range(3):
            if (x, y) != (1, 1):
                faces.extend(
                    [
                        [[x, y], [x + 1, y], [x, y + 1]],
                        [[x + 1, y], [x + 1, y + 1], [x, y + 1]],
                    ]
                )
    edges = boundary_segments(np.array(faces, dtype=float))
    assert len(edges) == 16
    assert sum(bool(((edge >= 1) & (edge <= 2)).all()) for edge in edges) == 4


def test_transparent_edges_do_not_fade_scalar_colors() -> None:
    """Change stroke opacity while retaining fully opaque scientific colormaps."""
    art = figure((4, 4))
    try:
        samples = _square()
        layer = LAYERS[2]
        collection = panel(art.subplots(), samples, layer, samples.values[0, 0], (0, 1))
        art.canvas.draw()
        np.testing.assert_allclose(np.asarray(collection.get_facecolor())[:, 3], 1)
        np.testing.assert_allclose(np.asarray(collection.get_edgecolor())[:, 3], 0.18)
        assert collection.cmap.name == "hot"
    finally:
        art.clear()


def test_convergence_fields_and_mesh_only_exports_are_separate(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Draw domain outlines on fields and the full neutral grid in a separate figure."""
    samples = _square()
    layer = replace(LAYERS[2], show_edges=False, outline=True)
    source = _source(tmp_path, convergence=True)
    inspected = []

    def inspect(art: Figure, destination: Path, options: RenderOptions) -> None:
        inspected.append(destination.name)
        collections = art.axes[0].collections
        if destination.name == "timecourse":
            assert len(collections) == 2
            assert not np.asarray(collections[0].get_edgecolor()).size
            outline = collections[1]
            assert isinstance(outline, LineCollection)
            assert len(outline.get_segments()) == 4
        else:
            assert len(collections) == 1
            assert collections[0].get_array() is None
            np.testing.assert_allclose(
                np.asarray(collections[0].get_facecolor()), [[1, 1, 1, 1]]
            )
            np.testing.assert_allclose(
                np.asarray(collections[0].get_edgecolor())[:, 3], 0.45
            )

    monkeypatch.setattr("research.spt.plot_render.save_figure", inspect)
    options = RenderOptions(layers=(layer,), snapshots=(0,))
    render_timecourse(source, samples, tmp_path, options, {layer.sid: (0, 1)})
    render_mesh(source, samples, tmp_path, options)
    assert inspected == ["timecourse", "mesh"]


@pytest.mark.parametrize("alpha", [-0.01, 1.01, float("nan")])
def test_invalid_edge_opacity_is_rejected(alpha: float) -> None:
    """Reject unsupported alpha values before producing partial figures."""
    with pytest.raises(ValueError, match="edge opacity"):
        RenderOptions(layers=(replace(LAYERS[0], edge_alpha=alpha),)).validate()
