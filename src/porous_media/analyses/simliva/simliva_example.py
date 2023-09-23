"""Example for plotting."""

from pathlib import Path
from typing import Dict, List

from porous_media import BASE_DIR
from porous_media.analyses.spt.spt_substrate_scan import visualize_scan
from porous_media.console import console
from porous_media.data.xdmf_tools import AttributeType, XDMFInfo, xdmfs_from_directory
from porous_media.log import get_logger
from porous_media.visualization.pyvista_visualization import DataLayer


logger = get_logger(__name__)

# Definition of variables for which panels are generated (this should be pretty
# complete and be created based on the XMDFInfo
data_layers_simliva: List[DataLayer] = [
    DataLayer(
        sid="rr_necrosis",
        title="Necrosis (0: alive, 1: death)",
        colormap="binary",
    ),
    DataLayer(
        sid="rr_(glc)",
        title="Glucose [mM]",
        color_limits=(2, 5),
    ),
    DataLayer(
        sid="rr_(o2)",
        title="Oxygen [mM]",
    ),
    DataLayer(
        sid="rr_(lac)",
        title="Lactate [mM]",
    ),
    DataLayer(
        sid="rr_(atp)",
        title="ATP [mM]",
    ),
    DataLayer(
        sid="rr_(nadh)",
        title="NADH [mM]",
    ),
    DataLayer(
        sid="pressure",
        title="Pressure [?]",
    ),
    DataLayer(
        sid="fluid_flux_TPM",
        title="Fluid flow [?/s]",
        colormap="RdBu",
        viz_type=AttributeType.VECTOR,
    ),
]

# subset of scalars to visualize
selection_simliva: List[str] = [
    "rr_necrosis",
    "rr_(glc)",
    "rr_(o2)",
    "rr_(lac)",
    "rr_(atp)",
    "rr_(nadh)",
    # "pressure",
]


if __name__ == "__main__":
    # FIXME: add additional variables to the VTKs
    # FIXME: position; external concentration; rates; ...

    # [1] process files
    # Put all simulations in a single folder and call the xdmfs_from_directory on the
    # folder
    xdmf_dir = Path("/home/mkoenig/git/porous_media/data/simliva_2023-09-20/")
    xdmf_dict: Dict[Path, Path] = xdmfs_from_directory(
        input_dir=Path("/home/mkoenig/git/porous_media/data/simliva_2023-09-20/"),
        xdmf_dir=xdmf_dir,
        overwrite=False,
    )
    xdmf_paths = list(xdmf_dict.keys())
    info: XDMFInfo = XDMFInfo.from_path(xdmf_path=xdmf_paths[0])
    console.print(info)

    # [2] Based on the variables in the dataset the DataLayers have to be created.
    # This is the subset of data for which visualizations (panels are generated)
    console.print(data_layers_simliva)

    # [3] create visualizations
    results_dir: Path = BASE_DIR / "results" / "simliva_2023-09-20"
    visualize_scan(
        xdmf_paths=xdmf_paths,
        data_layers=data_layers_simliva,
        results_dir=results_dir,
        selection=selection_simliva,
        create_panels=True,
    )
