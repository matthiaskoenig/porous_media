"""Compose annotated manuscript map panels using canonical SPT scalar styles."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import numpy as np
from matplotlib.axes import Axes
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.cm import ScalarMappable
from matplotlib.figure import Figure, SubFigure

from research.spt.fields import LAYERS, Layer
from research.spt.plot_data import PlotSource
from research.spt.plot_render import (
    RenderOptions,
    Samples,
    display_limits,
    figure,
    panel,
    scalar_bar,
)
from research.spt.publication_data import (
    PATTERN_NAMES,
    SIMULATION_TIMES,
    PublicationData,
)

FIELD_TITLES = ("Substrate", "Product", "Protein", "Toxic\ncompound", "Necrosis")


def save_publication(art: Figure, destination: Path, options: RenderOptions) -> None:
    """Keep annotations as vector text and bound dense map exports at the chosen DPI."""
    for extension in options.formats:
        art.savefig(
            destination.with_suffix(f".{extension}"),
            dpi=options.dpi,
            bbox_inches="tight",
        )


def map_panel(
    parent: Figure | SubFigure,
    axis: Axes,
    samples: Samples,
    layer: Layer,
    values: np.ndarray,
    ranges: dict[str, tuple[float, float]],
    scale: float,
    dpi: int,
) -> None:
    """Use the existing maps and five-label bars, scaling type with panel size."""
    collection = raster_panel(
        axis,
        samples,
        replace(layer, line_width=layer.line_width * scale),
        values,
        display_limits(layer, ranges),
        side=4 * scale,
        dpi=dpi,
    )
    scalar_bar(parent, axis, collection, layer)
    bar = parent.axes[-1]
    bar.set_title(layer.title, fontsize=12 * scale, fontweight="bold", pad=8 * scale)
    bar.tick_params(labelsize=10 * scale, length=2 * scale, pad=2 * scale)
    bar.xaxis.get_offset_text().set_fontsize(10 * scale)


def raster_panel(
    axis: Axes,
    samples: Samples,
    layer: Layer,
    values: np.ndarray,
    limits: tuple[float, float],
    *,
    side: float,
    dpi: int,
    mesh: bool = False,
) -> ScalarMappable:
    """Render one bounded tile, avoiding full-figure raster buffers per PDF panel."""
    tile = Figure(figsize=(side, side), dpi=dpi)
    canvas = FigureCanvasAgg(tile)
    try:
        tile_axis = tile.add_axes((0, 0, 1, 1))
        collection = panel(tile_axis, samples, layer, values, limits, mesh=mesh)
        mappable = ScalarMappable(norm=collection.norm, cmap=collection.cmap)
        canvas.draw()
        pixels = np.asarray(canvas.buffer_rgba()).copy()
        extent = (*tile_axis.get_xlim(), *tile_axis.get_ylim())
    finally:
        tile.clear()
    axis.imshow(pixels, extent=extent, origin="upper", interpolation="none")
    axis.set_aspect("equal")
    axis.set_axis_off()
    return mappable


def time_grid(
    parent: Figure | SubFigure,
    samples: Samples,
    ranges: dict[str, tuple[float, float]],
    *,
    scale: float,
    dpi: int,
) -> None:
    """Restore field headers, exact snapshot labels and the original downward arrow."""
    axes = parent.subplots(len(samples.times), len(LAYERS), squeeze=False)
    for row, time in enumerate(samples.times):
        values = samples.frame(float(time))
        for column, layer in enumerate(LAYERS):
            axis = axes[row, column]
            map_panel(parent, axis, samples, layer, values[column], ranges, scale, dpi)
            if row == 0:
                axis.set_title(
                    FIELD_TITLES[column],
                    fontsize=18 * scale,
                    fontweight="bold",
                    pad=12 * scale,
                )
            if column == 0:
                axis.text(
                    -0.05,
                    0.5,
                    f"{time / 3600:g} h",
                    transform=axis.transAxes,
                    ha="right",
                    va="center",
                    fontsize=16 * scale,
                    fontweight="bold",
                )
    middle = axes[len(samples.times) // 2, 0]
    middle.annotate(
        "",
        xy=(-0.42, -1),
        xytext=(-0.42, 1),
        xycoords="axes fraction",
        arrowprops={"arrowstyle": "-|>", "color": "black", "lw": 2 * scale},
        annotation_clip=False,
    )
    middle.text(
        -0.53,
        0.0,
        "Time",
        transform=middle.transAxes,
        rotation=90,
        ha="center",
        va="center",
        fontsize=16 * scale,
        fontweight="bold",
    )


def render_example(
    data: PublicationData,
    source: PlotSource,
    ranges: dict[str, tuple[float, float]],
    directory: Path,
    options: RenderOptions,
) -> None:
    """Update Figure 4's representative case at 0, 6, 12, 18 and 24 hours."""
    art = figure((15, 15))
    try:
        time_grid(
            art,
            data.snapshots(source, SIMULATION_TIMES),
            ranges,
            scale=0.7,
            dpi=options.dpi,
        )
        save_publication(art, directory / "figure-04", options)
    finally:
        art.clear()


def render_patterns(
    data: PublicationData,
    sources: list[PlotSource],
    ranges: dict[str, tuple[float, float]],
    directory: Path,
    options: RenderOptions,
) -> None:
    """Update Figure 5's two-column, three-row comparison with A-F annotations."""
    art = figure((18, 27))
    try:
        parents = art.subfigures(3, 2, wspace=0.01, hspace=0.02)
        assert isinstance(parents, np.ndarray)
        for pattern, (parent, source) in enumerate(
            zip(parents.flat, sources, strict=True)
        ):
            parent.suptitle(
                f"{chr(65 + pattern)}  {PATTERN_NAMES[pattern]}",
                fontsize=18,
                fontweight="bold",
            )
            time_grid(
                parent,
                data.snapshots(source, SIMULATION_TIMES),
                ranges,
                scale=0.37,
                dpi=options.dpi,
            )
        save_publication(art, directory / "figure-05", options)
    finally:
        art.clear()


def render_necrosis(
    data: PublicationData,
    sources: list[PlotSource],
    directory: Path,
    options: RenderOptions,
) -> None:
    """Update Figure 7 with named rows, ordered doses and the inflow direction."""
    art = figure((18, 14))
    outer = art.add_gridspec(2, 1, height_ratios=(0.45, 12))
    heading = art.add_subplot(outer[0])
    heading.set_axis_off()
    heading.text(
        0.5,
        0.8,
        "Substrate / drug inflow",
        ha="center",
        va="center",
        fontsize=16,
        fontweight="bold",
        transform=heading.transAxes,
    )
    heading.annotate(
        "",
        xy=(0.8, 0.25),
        xytext=(0.2, 0.25),
        xycoords="axes fraction",
        arrowprops={"arrowstyle": "-|>", "lw": 2, "color": "black"},
    )
    grid = outer[1].subgridspec(6, 8)
    layer = replace(LAYERS[-1], limits=(0, 1))
    try:
        for source in sources:
            condition = source.run.condition
            assert condition is not None
            row, column = condition.pattern_key, condition.boundary_flow_key
            axis = art.add_subplot(grid[row, column])
            samples = data.snapshots(source, (86400.0,))
            map_panel(
                art,
                axis,
                samples,
                layer,
                samples.values[0, -1],
                {layer.sid: (0, 1)},
                0.45,
                options.dpi,
            )
            if row == 0:
                axis.set_title(f"Dose {column + 1}", fontsize=12, fontweight="bold")
            if column == 0:
                axis.text(
                    -0.1,
                    0.5,
                    f"{chr(65 + row)}  {PATTERN_NAMES[row].replace(' ', chr(10), 1)}",
                    transform=axis.transAxes,
                    rotation=90,
                    ha="right",
                    va="center",
                    fontsize=12,
                    fontweight="bold",
                )
        save_publication(art, directory / "figure-07", options)
    finally:
        art.clear()
