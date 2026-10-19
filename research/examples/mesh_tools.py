"""Archived example; configure research paths before running."""

from pathlib import Path

from porous_media.console import console
from porous_media.data.xdmf_tools import XDMFInfo

if __name__ == "__main__":
    # simple example for adding variables (load interpolated xdmf)
    from research import RESULTS_DIR

    xdmf_path: Path = (
        RESULTS_DIR
        / "spt_zonation_patterns_219/10_28800.0/simulation_pattern1_interpolated.xdmf"
    )
    xdmf_info = XDMFInfo.from_path(xdmf_path)
    console.print(xdmf_info)

    # console.rule(style="white")
    # xr_cells: xr.Dataset
    # xr_points: xr.Dataset
    # xr_cells, xr_points = create_mesh_dataframes(xdmf_path)
    # console.print(xr_cells)
