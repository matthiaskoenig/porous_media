"""Adapt the canonical manuscript layers for historical PyVista workflows."""

from typing import cast

from pyvista.plotting._typing import ColormapOptions

from porous_media.visualization.pyvista_visualization import DataLayer
from research.spt.fields import SPT_LAYERS

data_layers_spt = [
    DataLayer(
        sid=layer.sid,
        title=layer.title,
        colormap=cast(ColormapOptions, layer.colormap),
        color_limits=layer.limits,
        conversion_factor=None if layer.factor == 1.0 else layer.factor,
        scalar_bar=layer.scalar_bar,
    )
    for layer in SPT_LAYERS
]
