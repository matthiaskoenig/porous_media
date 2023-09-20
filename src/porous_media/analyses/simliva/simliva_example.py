"""Example for plotting."""

from pathlib import Path
from typing import List, Dict

from porous_media import BASE_DIR
from porous_media.analyses.spt.spt_substrate_scan import visualize_scan
from porous_media.console import console
from porous_media.data.xdmf_tools import XDMFInfo, xdmfs_from_directory
from porous_media.log import get_logger
from porous_media.visualization.pyvista_visualization import DataLayer

logger = get_logger(__name__)

data_layers_simliva: List[DataLayer] = [
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

# subset of scalars to visualize
selection_simliva: List[str] = [
    "rr_protein",
    "rr_(S)",
    "rr_(P)",
    "rr_(T)",
    "rr_necrosis",
    # "pressure",
]


if __name__ == "__main__":
    # process files
    xdmf_dir = Path("/home/mkoenig/git/porous_media/data/simliva_2023-09-20/")
    xdmf_dict: Dict[Path, Path] = xdmfs_from_directory(
        input_dir=Path(
            "/home/mkoenig/git/porous_media/data/simliva_2023-09-20/"
        ),
        xdmf_dir=xdmf_dir,
        overwrite=False,
    )
    xdmf_paths = list(xdmf_dict.keys())
    info: XDMFInfo = XDMFInfo.from_path(xdmf_path=xdmf_paths[0])
    console.print(info)

    # create visualizations
    # from porous_media.analyses.spt import data_layers_spt
    # results_dir: Path = BASE_DIR / "results" / "simliva_2023-09-20"
    # visualize_scan(
    #     xdmf_paths=xdmf_paths,
    #     data_layers=data_layers_spt,
    #     results_dir=results_dir,
    #     selection=selection_simliva,
    #     create_panels=True,
    # )
