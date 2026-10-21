"""Scale manuscript exports while preserving their composition and export DPI."""

from __future__ import annotations

from typing import cast

import numpy as np
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.collections import Collection
from matplotlib.figure import Figure, SubFigure
from matplotlib.layout_engine import ConstrainedLayoutEngine
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, Patch
from matplotlib.text import Annotation, Text
from matplotlib.transforms import Bbox
from PIL import Image

from research.spt.publication_patterns import PatternIcon


def publication_bounds(art: Figure) -> Bbox:
    """Crop to displayed panels and annotations, with a one-point safety border.

    Hidden map axes must not contribute their empty rectangular slots. Include
    the actual image bounds, all annotations, subfigure titles and legends.
    Ordinary chart axes retain their ticks, labels and full plotting area.
    """
    art.canvas.draw()
    renderer = cast(FigureCanvasAgg, art.canvas).get_renderer()
    # Outside legends otherwise re-anchor themselves when savefig crops the
    # figure bounding box. Keep their measured position in figure coordinates.
    for legend in art.legends:
        box = legend.get_window_extent(renderer)
        x, y = art.transFigure.inverted().transform((box.x0, box.y0))
        legend.set_bbox_to_anchor((float(x), float(y)), transform=art.transFigure)
        legend.set_loc("lower left")
        legend.borderaxespad = 0
    bounds = []
    for axis in art.axes:
        if not axis.get_visible():
            continue
        if axis.axison:
            box = axis.get_tightbbox(renderer)
            if box is not None:
                bounds.append(box)
        else:
            for artist in (*axis.images, *axis.texts, axis.title):
                if artist.get_visible() and (
                    not isinstance(artist, Text)
                    or artist.get_text()
                    or isinstance(artist, Annotation)
                ):
                    bounds.append(artist.get_window_extent(renderer))
    parents: list[Figure | SubFigure] = [art]
    for parent in parents:
        parents.extend(parent.subfigs)
        for artist in parent.texts:
            if artist.get_visible():
                bounds.append(artist.get_window_extent(renderer))
        for legend in parent.legends:
            if not legend.get_visible():
                continue
            artists = (
                [legend]
                if legend.get_frame_on()
                else [*legend.get_texts(), legend.get_title(), *legend.legend_handles]
            )
            for artist in artists:
                if artist is not None and artist.get_visible():
                    bounds.append(artist.get_window_extent(renderer))
    if not bounds:
        raise ValueError("Publication figure contains no visible content")
    bounds.extend(
        icon.get_window_extent(renderer)
        for icon in art.findobj(match=PatternIcon)
        if icon.get_visible()
    )
    # Freeze the computed positions across backends so vector font metrics cannot
    # move panels outside the crop calculated by the Agg renderer.
    art.set_layout_engine("none")
    return Bbox.union(bounds).transformed(art.dpi_scale_trans.inverted()).padded(1 / 72)


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
        if isinstance(artist, PatternIcon):
            artist.scale(scale)
        elif isinstance(artist, Text):
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
    layout = art.get_layout_engine()
    if isinstance(layout, ConstrainedLayoutEngine):
        settings = layout.get()
        if "w_pad" in settings and "h_pad" in settings:
            layout.set(w_pad=settings["w_pad"] * scale, h_pad=settings["h_pad"] * scale)
