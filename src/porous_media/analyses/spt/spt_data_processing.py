"""Process all SPT simulations."""

from pathlib import Path
from typing import Dict, List

from porous_media.console import console
from porous_media.data.xdmf_calculations import mesh_datasets_from_xdmf
from porous_media.data.xdmf_tools import XDMFInfo, xdmfs_from_directory


def process_spt_simulations(input_dir: Path, xdmf_dir: Path) -> Dict[Path, Path]:
    """Process SPT simulation results."""
    # process files

    xdmfs: Dict[Path, Path] = xdmfs_from_directory(
        input_dir=input_dir,
        xdmf_dir=xdmf_dir,
        overwrite=False,
    )
    return xdmfs


def process_spt_meshes(xdmf_dir: Path) -> None:
    """Process the meshes for the SPT."""
    # prepare datasets for analysis
    xdmf_paths = sorted([p for p in xdmf_dir.glob("*.xdmf")])
    console.print(xdmf_paths)
    for xdmf_path in xdmf_paths:
        _, _ = mesh_datasets_from_xdmf(xdmf_path, overwrite=False)


if __name__ == "__main__":
    input_dir = Path("/data/qualiperf/P7-Perf/spt_results/simulation")
    xdmf_dir = Path("/home/mkoenig/git/porous_media/data/spt/2023-12-13/xdmf")

    # process_spt_simulations(input_dir=input_dir, output_dir=output_dir)
    process_spt_meshes(xdmf_dir=xdmf_dir)
