"""Verify GeoJSON plotting works without private data or optional packages."""

import json
from pathlib import Path

import matplotlib.pyplot as plt
import pytest

from porous_media.mesh.mesh_geojson import plot_geojson


def test_plot_geojson(tmp_path: Path) -> None:
    """Render a polygon and a multipolygon from an explicit local input."""
    ring = [[0, 0], [2, 0], [2, 2], [0, 0]]
    path = tmp_path / "lobuli.geojson"
    path.write_text(
        json.dumps(
            {
                "type": "FeatureCollection",
                "features": [
                    {
                        "type": "Feature",
                        "geometry": {"type": "Polygon", "coordinates": [ring]},
                    },
                    {
                        "type": "Feature",
                        "geometry": {"type": "MultiPolygon", "coordinates": [[ring]]},
                    },
                    {"type": "Feature", "geometry": None},
                ],
            }
        ),
        encoding="utf-8",
    )
    figure = plot_geojson(path)
    try:
        assert len(figure.axes[0].patches) == 2
        figure.savefig(tmp_path / "lobuli.png")
        assert (tmp_path / "lobuli.png").stat().st_size > 0
    finally:
        plt.close(figure)


def test_plot_geojson_requires_feature_collection(tmp_path: Path) -> None:
    """Reject an unsupported document before opening a figure."""
    path = tmp_path / "point.geojson"
    path.write_text('{"type": "Point", "coordinates": [0, 0]}', encoding="utf-8")
    with pytest.raises(ValueError, match="FeatureCollection"):
        plot_geojson(path)
