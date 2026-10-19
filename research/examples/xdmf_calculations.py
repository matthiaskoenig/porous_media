"""Archived example; configure research paths before running."""

from pathlib import Path

from porous_media.console import console
from porous_media.data.xdmf_calculations import mesh_datasets_from_xdmf

if __name__ == "__main__":
    console.rule(title="XDMF calculations", style="white")

    # prepare data for analysis
    xdmf_path = (
        Path("/home/mkoenig/git/porous_media/data/spt/2023-12-05") / "sim001.xdmf"
    )
    xr_cells, xr_points = mesh_datasets_from_xdmf(xdmf_path, overwrite=False)

    console.rule(align="left", title="cell_data", style="white")
    console.print(xr_cells)
    console.rule(align="left", title="point_data", style="white")
    console.print(xr_points)
