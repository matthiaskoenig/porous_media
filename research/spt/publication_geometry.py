"""Compose annotated manuscript map panels using canonical SPT scalar styles."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import numpy as np
from matplotlib.axes import Axes
from matplotlib.backend_bases import RendererBase
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.cm import ScalarMappable
from matplotlib.figure import Figure, SubFigure
from matplotlib.patches import Patch
from matplotlib.transforms import Bbox

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
    SIMULATION_TIMES,
    PublicationData,
    pattern_label,
)
from research.spt.publication_layout import publication_bounds, resize_publication
from research.spt.publication_patterns import add_pattern_icon

FIELD_TITLES = ("Substrate", "Product", "Protein", "Toxic\ncompound", "Necrosis")


def label_pattern_row(axis: Axes, pattern: int, *, offset: float = 65) -> None:
    """Label one pattern row with the shared A-F name and typography."""
    label = axis.annotate(
        pattern_label(pattern),
        xy=(0, 0.5),
        xycoords="axes fraction",
        xytext=(-offset, 0),
        textcoords="offset points",
        rotation=90,
        ha="center",
        va="center",
        fontsize=12,
        fontweight="bold",
        annotation_clip=False,
    )
    add_pattern_icon(label, pattern)


def shared_colorbar(
    parent: Figure | SubFigure,
    axis: Axes,
    collection: ScalarMappable,
    layer: Layer,
    *,
    scale: float,
    vertical: bool = False,
) -> None:
    """Reuse canonical normalization and scalar labels in one shared field bar."""
    scalar_bar(parent, axis, collection, layer, cax=axis, vertical=vertical)
    original_locator = axis.get_axes_locator()

    def inset_bar(_target: object, renderer: RendererBase) -> Bbox:
        bounds = (
            original_locator(axis, renderer)
            if original_locator is not None
            else axis.get_position(original=True)
        )
        return (
            bounds.shrunk(1, 0.8).anchored("C", bounds)
            if vertical
            else bounds.shrunk(0.8, 1).anchored("C", bounds)
        )

    axis.set_axes_locator(inset_bar)
    if vertical:
        axis.set_ylabel(layer.title, fontsize=max(10, 12 * scale), fontweight="bold")
    else:
        axis.set_title(
            layer.title, fontsize=12 * scale, fontweight="bold", pad=6 * scale
        )
    axis.tick_params(
        labelsize=max(9, 10 * scale) if vertical else 10 * scale,
        length=2 * scale,
        pad=2 * scale,
    )
    dimension = axis.yaxis if vertical else axis.xaxis
    dimension.get_offset_text().set_fontsize(
        max(9, 10 * scale) if vertical else 10 * scale
    )
    for item in axis.collections:
        item.set_rasterized(False)
        item.set_edgecolor("face")


def save_publication(
    art: Figure,
    destination: Path,
    options: RenderOptions,
    *,
    width_mm: float | None = None,
) -> None:
    """Keep annotations as vector text and bound dense map exports at the chosen DPI."""
    resize_publication(art, width_mm)
    art.set_dpi(options.dpi)
    bounds = publication_bounds(art)
    for extension in options.formats:
        art.savefig(
            destination.with_suffix(f".{extension}"),
            dpi=options.dpi,
            bbox_inches=bounds,
            pad_inches=0,
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
    *,
    sectors: bool = False,
    colorbar: bool = True,
    tight: bool = False,
) -> ScalarMappable:
    """Draw a canonical map, optionally adding its own five-label scalar bar."""
    collection = raster_panel(
        axis,
        samples,
        replace(layer, line_width=layer.line_width * scale),
        values,
        display_limits(layer, ranges),
        side=4 * scale,
        dpi=dpi,
        sectors=sectors,
        tight=tight,
    )
    if not colorbar:
        return collection
    scalar_bar(parent, axis, collection, layer)
    bar = parent.axes[-1]
    bar.set_title(layer.title, fontsize=12 * scale, fontweight="bold", pad=8 * scale)
    bar.tick_params(labelsize=10 * scale, length=2 * scale, pad=2 * scale)
    bar.xaxis.get_offset_text().set_fontsize(10 * scale)
    return collection


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
    sectors: bool = False,
    tight: bool = False,
) -> ScalarMappable:
    """Render one bounded tile, avoiding full-figure raster buffers per PDF panel."""
    tile = Figure(figsize=(side, side), dpi=dpi)
    canvas = FigureCanvasAgg(tile)
    try:
        tile_axis = tile.add_axes((0, 0, 1, 1))
        collection = panel(
            tile_axis, samples, layer, values, limits, mesh=mesh, sectors=sectors
        )
        mappable = ScalarMappable(norm=collection.norm, cmap=collection.cmap)
        canvas.draw()
        # Crop the tile to the domain, including white living tissue. Its scalar
        # colors must not determine whether a scientific region is kept.
        low, high = samples.faces.min(axis=(0, 1)), samples.faces.max(axis=(0, 1))
        projected = tile_axis.transData.transform(np.stack((low, high)))
        buffer = np.asarray(canvas.buffer_rgba())
        height, width = buffer.shape[:2]
        padding = max(1, int(np.ceil(layer.line_width * dpi / 72)))
        left, bottom = np.maximum(0, np.floor(projected[0]).astype(int) - padding)
        right, top = np.minimum(
            (width, height), np.ceil(projected[1]).astype(int) + padding
        )
        pixels = buffer[height - top : height - bottom, left:right].copy()
        corners = tile_axis.transData.inverted().transform(
            ((left, bottom), (right, top))
        )
        extent = (corners[0, 0], corners[1, 0], corners[0, 1], corners[1, 1])
        xlim, ylim = tile_axis.get_xlim(), tile_axis.get_ylim()
    finally:
        tile.clear()
    axis.imshow(pixels, extent=extent, origin="upper", interpolation="none")
    axis.set_xlim(*(extent[:2] if tight else xlim))
    axis.set_ylim(*(extent[2:] if tight else ylim))
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
    colorbars: bool = True,
    layers: tuple[Layer, ...] = LAYERS,
) -> tuple[ScalarMappable, ...]:
    """Draw field rows and time columns with optional shared bars on the right."""
    grid = parent.add_gridspec(
        len(layers) + 1,
        len(samples.times) + int(colorbars),
        height_ratios=[0.13] + [1] * len(layers),
        width_ratios=[1] * len(samples.times) + ([0.045] if colorbars else []),
        hspace=0.02,
        wspace=0.02,
    )
    header = parent.add_subplot(grid[0, : len(samples.times)])
    header.set_axis_off()
    font_size = max(12, 16 * scale)
    header.text(
        0.02,
        0.5,
        "Time",
        transform=header.transAxes,
        ha="left",
        va="center",
        fontsize=font_size,
        fontweight="bold",
    )
    header.annotate(
        "",
        xy=(0.98, 0.5),
        xytext=(0.15, 0.5),
        xycoords="axes fraction",
        arrowprops={"arrowstyle": "-|>", "color": "black", "lw": 2 * scale},
    )
    collections = []
    frames = tuple(samples.frame(float(time)) for time in samples.times)
    for row, layer in enumerate(layers):
        field = LAYERS.index(layer)
        for column, (time, values) in enumerate(
            zip(samples.times, frames, strict=True)
        ):
            axis = parent.add_subplot(grid[row + 1, column])
            collection = map_panel(
                parent,
                axis,
                samples,
                layer,
                values[field],
                ranges,
                scale,
                dpi,
                colorbar=False,
                tight=True,
            )
            if column == 0:
                collections.append(collection)
                axis.annotate(
                    FIELD_TITLES[field].replace("\n", " "),
                    xy=(0, 0.5),
                    xycoords="axes fraction",
                    xytext=(-10, 0),
                    textcoords="offset points",
                    rotation=90,
                    ha="center",
                    va="center",
                    fontsize=font_size,
                    fontweight="bold",
                    annotation_clip=False,
                )
            if row == 0:
                axis.set_title(
                    f"{time / 3600:g} h",
                    fontsize=font_size,
                    fontweight="bold",
                    pad=6 * scale,
                )
        if colorbars:
            shared_colorbar(
                parent,
                parent.add_subplot(grid[row + 1, -1]),
                collections[-1],
                layer,
                scale=scale,
                vertical=True,
            )
    return tuple(collections)


def render_example(
    data: PublicationData,
    source: PlotSource,
    ranges: dict[str, tuple[float, float]],
    directory: Path,
    options: RenderOptions,
    *,
    width_mm: float | None = None,
) -> None:
    """Show the representative case at five true times with shared field colorbars."""
    art = figure((12, 13))
    try:
        assert source.run.condition is not None and source.dose is not None
        heading = art.suptitle(
            f"{pattern_label(source.run.condition.pattern_key)} | D{source.dose + 1}",
            fontsize=16,
            fontweight="bold",
        )
        add_pattern_icon(heading, source.run.condition.pattern_key, size=24)
        time_grid(
            art,
            data.snapshots(source, SIMULATION_TIMES),
            ranges,
            scale=0.7,
            dpi=options.dpi,
        )
        save_publication(art, directory / "figure-04", options, width_mm=width_mm)
    finally:
        art.clear()


def render_patterns(
    data: PublicationData,
    sources: list[PlotSource],
    ranges: dict[str, tuple[float, float]],
    directory: Path,
    options: RenderOptions,
    *,
    width_mm: float | None = None,
) -> None:
    """Compare protein, toxic compound and necrosis with time running left to right."""
    art = figure((14, 16))
    try:
        assert sources[0].dose is not None
        art.suptitle(f"D{sources[0].dose + 1}", fontsize=16, fontweight="bold")
        outer = art.add_gridspec(
            4, 6, height_ratios=(1, 1, 1, 0.025), hspace=0.025, wspace=0.025
        )
        collections = ()
        for pattern, source in enumerate(sources):
            column = 3 * (pattern % 2)
            parent = art.add_subfigure(outer[pattern // 2, column : column + 3])
            heading = parent.suptitle(
                pattern_label(pattern),
                fontsize=16,
                fontweight="bold",
            )
            add_pattern_icon(heading, pattern, size=24)
            collections = time_grid(
                parent,
                data.snapshots(source, SIMULATION_TIMES),
                ranges,
                scale=0.37,
                dpi=options.dpi,
                colorbars=False,
                layers=LAYERS[2:],
            )
        for column, (layer, collection) in enumerate(
            zip(LAYERS[2:], collections, strict=True)
        ):
            shared_colorbar(
                art,
                art.add_subplot(outer[-1, 2 * column : 2 * column + 2]),
                collection,
                layer,
                scale=1,
            )
        save_publication(art, directory / "figure-05", options, width_mm=width_mm)
    finally:
        art.clear()


def render_necrosis(
    data: PublicationData,
    sources: list[PlotSource],
    directory: Path,
    options: RenderOptions,
) -> None:
    """Compare A-F rows and D1-D8 columns with one shared binary state legend."""
    art = figure((16, 12))
    outer = art.add_gridspec(
        2, 2, height_ratios=(0.45, 12), width_ratios=(24, 0.6), wspace=0.02
    )
    heading = art.add_subplot(outer[0, 0])
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
    grid = outer[1, 0].subgridspec(6, 8, hspace=0.02, wspace=0.02)
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
                colorbar=False,
                tight=True,
            )
            if row == 0:
                axis.set_title(f"D{column + 1}", fontsize=12, fontweight="bold")
            if column == 0:
                label_pattern_row(axis, row, offset=14)
        legend_axis = art.add_subplot(outer[1, 1])
        legend_axis.set_axis_off()
        art.legend(
            handles=[
                Patch(facecolor="black", edgecolor="black", label="necrosis"),
                Patch(facecolor="white", edgecolor="black", label="no necrosis"),
            ],
            loc="center left",
            bbox_to_anchor=(0, 0.5),
            bbox_transform=legend_axis.transAxes,
            frameon=False,
            fontsize=14,
            handlelength=1.2,
            handleheight=1.2,
            labelspacing=1.1,
            borderaxespad=0,
        )
        save_publication(art, directory / "figure-07", options)
    finally:
        art.clear()
