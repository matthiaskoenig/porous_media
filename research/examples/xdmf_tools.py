"""Archived example; configure research paths before running."""

from pathlib import Path

from porous_media.console import console
from porous_media.data.xdmf_tools import DataLimits, XDMFInfo, vtks_to_xdmf

if __name__ == "__main__":
    from porous_media import RESOURCES_DIR
    from research import RESULTS_DIR

    # Simple test that all variables are read
    # vtk_dir = RESOURCES_DIR / "vtk" / "vtk_single"
    vtk_dir = RESOURCES_DIR / "vtk" / "vtk_timecourse"
    xdmf_path = RESULTS_DIR / "vtk_test.xdmf"
    vtks_to_xdmf(vtk_dir, xdmf_path=xdmf_path, overwrite=True)
    xdmf_info: XDMFInfo = XDMFInfo.from_path(xdmf_path)
    console.print(xdmf_info)

    xdmf_paths: list[Path] = [
        Path(f"/home/mkoenig/git/porous_media/data/spt_substrate_scan/sim_{k}.xdmf")
        for k in range(21, 26)
    ]

    xdmf_path = Path(
        "/home/mkoenig/git/porous_media/data/spt_substrate_scan/sim_25.xdmf"
    )
    xdmf_info = XDMFInfo.from_path(xdmf_path)
    console.print(xdmf_info)

    # limits of individual simulations
    all_limits: list[DataLimits] = []
    for xdmf_path in xdmf_paths:
        limits = DataLimits.from_xdmf(xdmf_path=xdmf_path, overwrite=False)
        all_limits.append(limits)

    # merge data limits
    data_limits = DataLimits.merge_limits(all_limits)
    console.print(data_limits)
