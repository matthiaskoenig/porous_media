"""Visualization-only SPT layer configuration, imported on demand."""

from porous_media.data.xdmf_tools import AttributeType
from porous_media.visualization.pyvista_visualization import DataLayer

# information for visualization
data_layers_spt: list[DataLayer] = [
    DataLayer(
        sid="rr_necrosis",
        title="Necrosis (0: alive, 1: death)",
        colormap="binary",
        viz_type=AttributeType.SCALAR,
    ),
    DataLayer(
        sid="rr_protein",
        title="Protein",
        colormap="hot",
    ),
    DataLayer(
        sid="rr_(S_ext)",
        title="Substrate S plasma [mM]",
        colormap="magma",
    ),
    DataLayer(
        sid="rr_(P_ext)",
        title="Product P plasma [mM]",
    ),
    DataLayer(
        sid="rr_(S)",
        title="Substrate S [mM]",
    ),
    DataLayer(
        sid="rr_(P)",
        title="Product P [mM]",
    ),
    DataLayer(
        sid="rr_(T)",
        title="Toxic compound T [mM]",
    ),
    DataLayer(
        sid="pressure",
        title="Pressure [?]",
    ),
]
