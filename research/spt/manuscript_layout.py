"""Scale manuscript exports while preserving their composition and export DPI."""

from __future__ import annotations

from typing import cast

from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure, SubFigure
from matplotlib.text import Annotation, Text
from matplotlib.transforms import Bbox

from research.spt.manuscript_patterns import PatternIcon


def manuscript_bounds(art: Figure) -> Bbox:
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
        raise ValueError("Manuscript figure contains no visible content")
    bounds.extend(
        icon.get_window_extent(renderer)
        for icon in art.findobj(match=PatternIcon)
        if icon.get_visible()
    )
    # Freeze the computed positions across backends so vector font metrics cannot
    # move panels outside the crop calculated by the Agg renderer.
    art.set_layout_engine("none")
    return Bbox.union(bounds).transformed(art.dpi_scale_trans.inverted()).padded(1 / 72)
