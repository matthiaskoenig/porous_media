"""Definition of SPT information."""

from pathlib import Path

from porous_media.data.xdmf_tools import AttributeType
from porous_media.visualization.pyvista_visualization import DataLayer


results_date: str = "2026-04-27"
data_spt_dir: Path = Path("/media/mkoenig/Extreme Pro/spt/data") / results_date
results_spt_dir: Path = Path("/media/mkoenig/Extreme Pro/spt/data") / results_date

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

# subset of scalars to visualize
selection_spt: list[str] = [
    "rr_(S_ext)",
    "rr_(P_ext)",
    "rr_protein",
    "rr_(T)",
    "rr_necrosis",
    # "pressure",
]
