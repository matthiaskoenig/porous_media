"""Zonation pattern scan."""

from pathlib import Path
from typing import List, Dict

from porous_media import BASE_DIR
from porous_media.analyses.spt.spt_substrate_scan import visualize_scan
from porous_media.console import console
from porous_media.data.xdmf_tools import XDMFInformation, xdmfs_from_directory
from porous_media.log import get_logger


logger = get_logger(__name__)

if __name__ == "__main__":
    # process files
    xdmf_dir = Path("/home/mkoenig/git/porous_media/data/spt/spt_zonation_patterns_new")
    xdmf_dict: Dict[Path, Path] = xdmfs_from_directory(
        input_dir=Path(
            "/home/mkoenig/git/porous_media/data/spt/simulation_zonation_hard_new"
        ),
        xdmf_dir=xdmf_dir,
        overwrite=False,
    )
    xdmf_paths = list(xdmf_dict.keys())
    info: XDMFInformation = XDMFInformation.from_path(xdmf_path=xdmf_paths[0])
    console.print(info)

    # create visualizations
    from porous_media.analyses.spt import data_layers_spt

    results_dir: Path = BASE_DIR / "results" / "spt_zonation_patterns_new"
    visualize_scan(
        xdmf_paths=xdmf_paths,
        data_layers=data_layers_spt,
        results_dir=results_dir,
        create_panels=True,
    )
