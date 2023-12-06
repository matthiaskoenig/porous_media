"""Process all SPT simulations."""

from pathlib import Path
from typing import Dict

from porous_media.console import console
from porous_media.data.xdmf_calculations import mesh_datasets_from_xdmf
from porous_media.data.xdmf_tools import (
    XDMFInfo,
    xdmfs_from_directory,
)


def process_spt_simulations(input_dir: Path, output_dir: Path) -> None:
    """Process SPT simulation results."""
    # process files

    xdmf_paths: Dict[Path, Path] = xdmfs_from_directory(
        input_dir=input_dir,
        xdmf_dir=output_dir,
        overwrite=False,
    )
    info: XDMFInfo = XDMFInfo.from_path(list(xdmf_paths.keys())[0])
    console.print(info)

    # prepare datasets for analysis
    for xdmf_path in xdmf_paths:
        _, _ = mesh_datasets_from_xdmf(xdmf_path, overwrite=False)


if __name__ == "__main__":

    input_dir = Path("/data/qualiperf/P7-Perf/spt_results/simulation_2711")
    output_dir = Path("/home/mkoenig/git/porous_media/data/spt/2023-12-05")
    process_spt_simulations(input_dir=input_dir, output_dir=output_dir)
