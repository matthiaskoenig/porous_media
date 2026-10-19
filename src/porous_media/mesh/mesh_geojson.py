"""Plot polygon boundaries from a GeoJSON feature collection.

Run ``python -m porous_media.mesh.mesh_geojson lobuli.geojson`` to inspect a
local file. Importing this module never opens files or starts a plot window.
"""

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.patches import PathPatch
from matplotlib.path import Path as PlotPath


def plot_geojson(geojson_path: Path) -> Figure:
    """Plot Polygon and MultiPolygon features, preserving interior rings.

    Args:
        geojson_path: Path to a GeoJSON FeatureCollection.

    Returns:
        Figure containing the polygon boundaries, with equal axis scaling.

    Raises:
        ValueError: If the document is not a FeatureCollection.
    """
    data = json.loads(geojson_path.read_text(encoding="utf-8"))
    if data.get("type") != "FeatureCollection":
        raise ValueError("Expected a GeoJSON FeatureCollection")

    figure, axes = plt.subplots(figsize=(10, 10))
    for feature in data["features"]:
        geometry = feature.get("geometry")
        if geometry is None:
            continue
        geometry_type = geometry["type"]
        if geometry_type == "Polygon":
            polygons = [geometry["coordinates"]]
        elif geometry_type == "MultiPolygon":
            polygons = geometry["coordinates"]
        else:
            continue
        for polygon in polygons:
            rings = [PlotPath([point[:2] for point in ring]) for ring in polygon]
            if rings:
                axes.add_patch(
                    PathPatch(
                        PlotPath.make_compound_path(*rings),
                        facecolor="#6699cc",
                        edgecolor="black",
                        alpha=0.5,
                    )
                )
    axes.autoscale_view()
    axes.set_aspect("equal")
    return figure


def main() -> None:
    """Plot a GeoJSON file supplied on the command line."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("geojson_path", type=Path)
    arguments = parser.parse_args()
    plot_geojson(arguments.geojson_path)
    plt.show()


if __name__ == "__main__":
    main()
