"""Tools to manipulate and work with meshes."""

from __future__ import annotations

from pathlib import Path

import meshio

from porous_media.console import console
from porous_media.data.xdmf_tools import XDMFInfo


def mesh_to_xdmf(m: meshio.Mesh, xdmf_path: Path, test_read: bool = False) -> None:
    """Serialize mesh to XDMF."""
    m.write(xdmf_path)
    if test_read:
        # check that the serialized mesh can be read again
        _ = meshio.read(xdmf_path)


def mesh_to_vtk(m: meshio.Mesh, vtk_path: Path, test_read: bool = False) -> None:
    """Serialize mesh to VTK.

    FIXME: must handle the special case of timecourse data; i.e. only serialze part of it
    """
    m.write(vtk_path)
    if test_read:
        # check that the serialized mesh can be read again
        _ = meshio.read(vtk_path)


if __name__ == "__main__":
    # simple example for adding variables (load interpolated xdmf)
    from porous_media import RESULTS_DIR

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
