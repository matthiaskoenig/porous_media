"""Tools to manipulate and work with meshes."""

from __future__ import annotations

from pathlib import Path

import meshio

from porous_media.console import console


def calculate_geometry_variables() -> None:
    """Calculate volumes and other important variables using meshplex.

    These will be added to the cell_data and point data.
    https://meshplex.readthedocs.io/
    """
    # FIXME: implement
    pass


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
