"""Scale manuscript exports while preserving their composition and export DPI."""

from __future__ import annotations

from typing import cast

import numpy as np
from matplotlib.collections import Collection
from matplotlib.figure import Figure
from matplotlib.layout_engine import ConstrainedLayoutEngine
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, Patch
from matplotlib.text import Annotation, Text
from PIL import Image


def resize_publication(art: Figure, width_mm: float | None) -> None:
    """Scale physical dimensions, typography and raster tiles together before export.

    Width describes the nominal canvas before tight cropping. Scalar arrays and
    line coordinates remain unchanged. Only embedded display images are resized.
    """
    if width_mm is None:
        return
    art.canvas.draw()
    scale = width_mm / (25.4 * art.get_figwidth())
    art.set_size_inches(*(art.get_size_inches() * scale))
    for artist in art.findobj():
        if isinstance(artist, Text):
            artist.set_fontsize(
                artist.get_fontproperties().get_size_in_points() * scale
            )
            if isinstance(artist, Annotation) and artist.anncoords == "offset points":
                x, y = artist.get_position()
                artist.set_position((x * scale, y * scale))
            if isinstance(artist, Annotation) and artist.arrow_patch is not None:
                arrow = artist.arrow_patch
                arrow.set_linewidth(arrow.get_linewidth() * scale)
                arrow.set_mutation_scale(arrow.get_mutation_scale() * scale)
        elif isinstance(artist, Line2D):
            artist.set_linewidth(artist.get_linewidth() * scale)
            artist.set_markersize(artist.get_markersize() * scale)
            artist.set_markeredgewidth(artist.get_markeredgewidth() * scale)
        elif isinstance(artist, Collection):
            artist.set_linewidth((np.asarray(artist.get_linewidth()) * scale).tolist())
        elif isinstance(artist, Patch):
            artist.set_linewidth(artist.get_linewidth() * scale)
            if isinstance(artist, FancyArrowPatch):
                artist.set_mutation_scale(artist.get_mutation_scale() * scale)
    for axis in art.axes:
        for dimension in (axis.xaxis, axis.yaxis):
            dimension.labelpad *= scale
            for which, ticks in (
                ("major", dimension.get_major_ticks()),
                ("minor", dimension.get_minor_ticks()),
            ):
                if ticks:
                    tick = ticks[0]
                    dimension.set_tick_params(
                        which=which,
                        length=tick.tick1line.get_markersize(),
                        width=tick.tick1line.get_markeredgewidth(),
                        pad=cast(float, tick.get_pad()) * scale,
                        labelsize=tick.label1.get_fontproperties().get_size_in_points(),
                    )
        for image in axis.images:
            pixels = image.get_array()
            if pixels is None or scale >= 1:
                continue
            with Image.fromarray(np.asarray(pixels)) as original:
                width, height = original.size
                size = max(1, round(width * scale)), max(1, round(height * scale))
                with original.resize(size, Image.Resampling.LANCZOS) as smaller:
                    image.set_data(np.asarray(smaller).copy())
    # The outside legends cache handle spacing in points when they are created.
    # Rebuild them at the new text size so that spacing scales with the canvas.
    for legend in list(art.legends):
        handles = legend.legend_handles
        labels = [text.get_text() for text in legend.get_texts()]
        title = legend.get_title()
        font = legend.get_texts()[0].get_fontproperties()
        legend.remove()
        art.legend(
            handles=handles,
            labels=labels,
            title=title.get_text(),
            title_fontproperties=title.get_fontproperties(),
            prop=font,
            loc="outside upper center",
            ncol=len(handles),
            frameon=False,
        )
    layout = art.get_layout_engine()
    if isinstance(layout, ConstrainedLayoutEngine):
        settings = layout.get()
        if "w_pad" in settings and "h_pad" in settings:
            layout.set(w_pad=settings["w_pad"] * scale, h_pad=settings["h_pad"] * scale)
