"""Tools to manipulate and work with meshes."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, Tuple

import meshio
import numpy as np

from porous_media.console import console


def mesh_to_xdmf(m: meshio.Mesh, xdmf_path: Path, test_read: bool = False) -> None:
    """Serialize mesh to XDMF."""
    console.rule(title="Mesh Serialization XDMF", style="white")
    console.print(f"{m=}")
    m.write(xdmf_path)
    if test_read:
        # check that the serialzed mesh can be read again
        _ = meshio.read(xdmf_path)


def mesh_to_vtk(m: meshio.Mesh, vtk_path: Path, test_read: bool = False) -> None:
    """Serialize mesh to VTK."""
    console.rule(title="Mesh Serialization VTK", style="white")
    console.print(f"{m=}")
    m.write(vtk_path)
    if test_read:
        # check that the serialzed mesh can be read again
        _ = meshio.read(vtk_path)


@dataclass
class MeshTimepoint:
    """Class for storing VTK simulation results."""

    mesh: meshio.Mesh
    cell_data: Dict[str, np.ndarray] = field(repr=False)
    cell_data_shape: Dict[str, np.ndarray]
    point_data: Dict[str, np.ndarray] = field(repr=False)
    point_data_shape: Dict[str, np.ndarray]

    @classmethod
    def from_vtk(cls, vtk_path: Path, show: bool = False) -> MeshTimepoint:
        """Read mesh information and data from single VTK.

        :param vtk_path: Path to VTK file.
        :param show: boolean flag to show information.
        """
        if show:
            console.rule(align="left", style="white")
            console.print(vtk_path)

        mesh: meshio.Mesh = meshio.read(vtk_path)

        cell_data: Dict[str, np.ndarray] = {}
        cell_data_shape: Dict[str, np.ndarray] = {}
        for key in mesh.cell_data.keys():
            data = cls.process_cell_data(key=key, mesh=mesh)
            cell_data[key] = data
            cell_data_shape[key] = data.shape

        point_data: Dict[str, np.ndarray] = {}
        point_data_shape: Dict[str, np.ndarray] = {}
        for key in mesh.point_data.keys():
            data = cls.process_point_data(key=key, mesh=mesh)
            point_data[key] = data
            point_data_shape[key] = data.shape

        mesh_results = MeshTimepoint(
            mesh=mesh,
            cell_data=cell_data,
            cell_data_shape=cell_data_shape,
            point_data=point_data,
            point_data_shape=point_data_shape,
        )
        if show:
            console.print(mesh_results)
        return mesh_results

    @staticmethod
    def process_cell_data(key: str, mesh: meshio.Mesh) -> np.ndarray:
        """Process data from the mesh."""
        data = mesh.cell_data[key]
        data = np.dstack(data)
        return data.squeeze()

    @staticmethod
    def process_point_data(key: str, mesh: meshio.Mesh) -> np.ndarray:
        """Process data from the mesh."""
        data = mesh.point_data[key]
        data = np.dstack(data)
        return data.squeeze()


if __name__ == "__main__":
    from porous_media import DATA_DIR, EXAMPLE_VTK

    vtk_path = EXAMPLE_VTK

    mesh_tp: MeshTimepoint = MeshTimepoint.from_vtk(vtk_path=vtk_path, show=True)
    console.print("cell_type")
    console.print(mesh_tp.cell_data["cell_type"])

    console.rule(style="white")

    # TODO: 1. plot mesh data for different timepoints (2d mesh)
    # 4 degree & 37 degree
    # glc, oxygen, lactate, atp, ros, cell death, alt, ast
    # t = 0, 30, 60, 90, 120, 150, 180

    # 'cell_type':
    #   0: internal node
    #   1: periportal (outside)
    #   2: perivenous (inside)
