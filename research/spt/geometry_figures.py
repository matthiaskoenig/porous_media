"""Render planar sector comparisons without a display server or VTK renderer."""

from __future__ import annotations

import tempfile
from pathlib import Path
from typing import cast

import matplotlib.pyplot as plt
import meshio
import numpy as np
from matplotlib.collections import PolyCollection
from matplotlib.colors import BoundaryNorm, ListedColormap

from research.spt.artifacts import fingerprint, generation_path, publish, run_lock
from research.spt.conditions import pattern_idx2name
from research.spt.fields import SPT_LAYERS
from research.spt.geometry import GeometryJob, geometry_manifest

COMPARISON_PLOT_VERSION = 2


def render_comparison(job: GeometryJob, descriptor: Path) -> None:
    """Save mesh/zonation and final-state panels beside a derived time series.

    These thin extruded meshes are viewed along z. One triangular face per wedge
    represents its cell value; separate sector nodes preserve discontinuities.
    """
    with meshio.xdmf.TimeSeriesReader(descriptor) as reader:
        points, cells = reader.read_points_cells()
        _, _, initial = reader.read_data(0)
        end, _, final = reader.read_data(reader.num_steps - 1)
    faces = []
    for block in cells:
        if block.type in {"wedge", "triangle"}:
            faces.extend(points[block.data[:, :3], :2])
        elif block.type == "quad":
            faces.extend(points[block.data, :2])
        else:
            raise ValueError(f"Unsupported planar figure cell type: {block.type}")
    mixed_patterns = job.study == "simulations_combinations"
    selections = (
        (
            ("rr_protein", "Initial protein zonation", initial, None),
            ("rr_necrosis", f"Necrosis at {end:g} s", final, None),
        )
        if mixed_patterns
        else (
            ("sector_id", "Convergence meshes", initial, "tab10"),
            ("rr_(T)", f"Toxicant at {end:g} s", final, None),
        )
    )
    figure, axes = plt.subplots(1, 2, figsize=(11, 6), layout="constrained")
    try:
        for axis, (field, title, data, color) in zip(axes, selections, strict=True):
            values = np.concatenate(data[field]).reshape(-1)
            if field == "sector_id":
                cmap = ListedColormap(plt.get_cmap(color)(np.arange(6)))
                collection = PolyCollection(
                    faces,
                    cmap=cmap,
                    norm=BoundaryNorm(np.arange(7) - 0.5, 6),
                    linewidths=0.15,
                    edgecolors="black",
                )
            else:
                layer = next(layer for layer in SPT_LAYERS if layer.sid == field)
                collection = PolyCollection(
                    faces,
                    cmap=layer.colormap,
                    linewidths=layer.line_width,
                    edgecolors=layer.edge_color,
                )
                if field == "rr_necrosis":
                    collection.set_clim(0, 1)
            collection.set_array(values)
            axis.add_collection(collection)
            extent = float(np.linalg.norm(points[:, :2], axis=1).max()) * 1.25
            axis.set_xlim(-extent, extent)
            axis.set_ylim(-extent, extent)
            axis.set_aspect("equal")
            axis.set_axis_off()
            axis.set_title(title)
            if field != "sector_id":
                figure.colorbar(collection, ax=axis, shrink=0.7)
            sector_ids = np.concatenate(initial["sector_id"]).reshape(-1)
            centers = np.concatenate(
                [points[block.data, :2].mean(axis=1) for block in cells]
            )
            for sector in range(6):
                position = centers[sector_ids == sector].mean(axis=0) * 1.8
                axis.text(
                    *position,
                    str(sector + 1),
                    ha="center",
                    va="center",
                    fontsize=11,
                    bbox={
                        "facecolor": "white",
                        "alpha": 0.85,
                        "edgecolor": "none",
                        "pad": 2,
                    },
                )
        if mixed_patterns:
            condition = job.sources[0].condition
            if condition is None:
                raise ValueError("Pattern comparison requires dose conditions")
            title = f"Dose {condition.boundary_flow_key + 1}: boundary flow {condition.boundary_flow:g} (source units unspecified)"
            labels = [
                f"{index + 1}: {pattern_idx2name[source.condition.pattern_key]}"
                for index, source in enumerate(job.sources)
                if source.condition
            ]
        else:
            title = "Six convergence resolutions: independent sectors"
            labels = [
                f"{index + 1}: {source.run_id.removeprefix('convergence_')}"
                for index, source in enumerate(job.sources)
            ]
        figure.suptitle(title)
        figure.supxlabel(" | ".join(labels), fontsize=9)
        figure.savefig(descriptor.with_suffix(".png"), dpi=200)
        figure.savefig(descriptor.with_suffix(".pdf"))
    finally:
        plt.close(figure)


def refresh_comparison(job: GeometryJob, output_root: Path) -> dict[str, object]:
    """Republish preview styling while retaining the existing simulation arrays.

    Hard links keep the large, immutable HDF sidecar unchanged. Publication uses
    the same atomic generation switch as reconstruction and preserves the previous
    generation if rendering fails.
    """
    directory = job.directory(output_root)
    with run_lock(directory):
        current = geometry_manifest(directory)
        if current.get("metadata", {}) == {
            "comparison_plot_version": COMPARISON_PLOT_VERSION
        }:
            return current
        previous = generation_path(directory, current)
        paths = tuple(
            Path(str(record["path"]))
            for record in cast(list[dict[str, object]], current["sources"])
        )
        before = [fingerprint(path) for path in paths]
        with tempfile.TemporaryDirectory(
            prefix=".partial-", dir=directory
        ) as temporary:
            staging = Path(temporary)
            for path in previous.iterdir():
                if path.suffix not in {".png", ".pdf"}:
                    (staging / path.name).hardlink_to(path)
            render_comparison(job, staging / f"{job.run_id}.xdmf")
            catalog = cast(dict[str, object], current["catalog"])
            publish(
                directory,
                staging,
                job.configuration(),
                paths,
                before,
                {"comparison_plot_version": COMPARISON_PLOT_VERSION},
                catalog=lambda generation: {
                    **catalog,
                    "xdmf": str(generation / f"{job.run_id}.xdmf"),
                },
            )
        return geometry_manifest(directory)
