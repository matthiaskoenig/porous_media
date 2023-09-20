"""Definition of SPT information."""
from typing import List

from porous_media.visualization.pyvista_visualization import DataLayer

data_layers_spt: List[DataLayer] = [
    DataLayer(
        sid="rr_necrosis",
        title="Necrosis (0: alive, 1: death)",
        colormap="binary",
        data_type="Scalar",
    ),
    DataLayer(
        sid="rr_protein",
        title="Protein",
        colormap="hot",
        data_type="Scalar",
    ),
    DataLayer(
        sid="rr_(S_ext)",
        title="Substrate S plasma [mM]",
        colormap="magma",
        data_type="Scalar",
    ),
    DataLayer(
        sid="rr_(P_ext)",
        title="Product P plasma [mM]",
        colormap="magma",
        data_type="Scalar",
    ),
    DataLayer(
        sid="rr_(S)",
        title="Substrate S [mM]",
        colormap="magma",
        data_type="Scalar",
    ),
    DataLayer(
        sid="rr_(P)",
        title="Product P [mM]",
        colormap="magma",
        data_type="Scalar",
    ),
    DataLayer(
        sid="rr_(T)",
        title="Toxic compound T [mM]",
        colormap="magma",
        data_type="Scalar",
    ),
    DataLayer(
        sid="pressure",
        title="Pressure [?]",
        colormap="magma",
        data_type="Scalar",
    ),
    # FIXME: visualize vector field
    DataLayer(
        sid="fluid_flux_TPM",
        title="Fluid flow [?/s]",
        colormap="RdBu",
        data_type="Vector",
    ),
]
