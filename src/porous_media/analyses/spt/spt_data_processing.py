"""Process all SPT simulations."""

from pathlib import Path
from typing import Dict, Iterable, List, Set

import numpy as np


from porous_media.console import console
from porous_media.data.xdmf_calculations import mesh_datasets_from_xdmf
from porous_media.data.xdmf_tools import (
    DataLimits,
    XDMFInfo,
    interpolate_xdmf,
    xdmfs_from_directory,
)
from porous_media.log import get_logger


logger = get_logger(__name__)


if __name__ == "__main__":

    # process files
    xdmf_dir = Path("/home/mkoenig/git/porous_media/data/spt/2023-12-05")
    xdmf_paths: Dict[Path, Path] = xdmfs_from_directory(
        input_dir=Path("/data/qualiperf/P7-Perf/spt_results/simulation"),
        xdmf_dir=xdmf_dir,
        overwrite=False,
    )
    info: XDMFInfo = XDMFInfo.from_path(list(xdmf_paths.keys())[0])
    console.print(info)

    # prepare datasets for analysis
    for xdmf_path in xdmf_paths:
        _, _ = mesh_datasets_from_xdmf(xdmf_path, overwrite=False)


    # interpolate data
    # Calculate tend time from all simulations
    # tends: np.ndarray = np.zeros(shape=(len(list(xdmf_paths)),))
    # for k, xdmf_path in enumerate(xdmf_paths):
    #     xdmf_info = XDMFInfo.from_path(xdmf_path)
    #     tends[k] = xdmf_info.tend
    # tend = tends.min()
    # for num in [10, 100]:
    #     output_dir = xdmf_dir / f"{num}"
    #     output_dir.mkdir(exist_ok=True, parents=True)
    #
    #     all_limits: List[DataLimits] = []
    #     for xdmf_path in xdmf_paths:
    #         interpolate_xdmf(
    #             xdmf_in=xdmf_path,
    #             xdmf_out=output_dir / f"{xdmf_path.stem}_interpolated.xdmf",
    #             times_interpolate=np.linspace(0, tend, num=num),
    #             overwrite=False,
    #         )
