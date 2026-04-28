"""Process all SPT simulations.

Results have to be processed into XDMF.
"""

from pathlib import Path

from porous_media.data.xdmf_tools import xdmfs_from_directory
from porous_media.data.xdmf_calculations import mesh_datasets_from_xdmf_dir


def process_febio_simulations(input_dir: Path, xdmf_dir: Path) -> dict[Path, Path]:
    """Process SPT simulation results."""
    # process files

    xdmfs: dict[Path, Path] = xdmfs_from_directory(
        input_dir=input_dir,
        xdmf_dir=xdmf_dir,
        overwrite=False,
    )
    return xdmfs


if __name__ == "__main__":
    from porous_media.analyses.spt import (
        data_spt_dir,
    )

    # [1] Process FEBIO simulations
    # for dir_name in ["convergence_sixth", "simulations_sixth"]:
    #     # process_spt_simulations(input_dir=input_dir, output_dir=output_dir)

    # [2] Hexagon reconstruction (parallel)
    # for in_dir, out_dir in [
    #     ("convergence_sixth", "convergence_lobulus"),
    #     # ("simulations_sixth", "simulations_lobulus"),
    # ]:
    #     # reconstruct lobulus from sixth
    #     xdmf_in_dir: Path = data_spt_dir / in_dir / "xdmf"
    #     xdmf_out_dir: Path = data_spt_dir / out_dir / "xdmf"
    #     reconstruct_lobulus_from_hexagon_dir(
    #         xdmf_in_dir=xdmf_in_dir,
    #         xdmf_out_dir=xdmf_out_dir,
    #     )

    # [3] Create mesh datasets
    for dir_name in [
        # "convergence_sixth",
        # "simulations_sixth",
        "convergence_lobulus",
        # "simulations_lobulus",
    ]:
        # process_spt_simulations(input_dir=input_dir, output_dir=output_dir)

        mesh_datasets_from_xdmf_dir(
            xdmf_dir=data_spt_dir / dir_name / "xdmf", overwrite=True
        )
