from __future__ import annotations

from pathlib import Path
from typing import Dict, Tuple, Iterable

import numpy as np
from pm.console import console
import meshio
from dataclasses import dataclass, field


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
        """Reads mesh information and data from single VTK.

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


@dataclass
class MeshTimecourse(MeshTimepoint):
    """Class for storing VTK simulation results."""

    time: np.ndarray
    time_shape: np.ndarray

    @classmethod
    def from_vtk_directory(cls, vtk_dir: Path, show: bool = False) -> MeshTimecourse:
        """Read MeshTimecourse from directory of VTKs."""
        if not vtk_dir.exists():
            raise IOError(f"VTK directory does not exist '{vtk_dir}'")

        console.print("vtk_paths:")
        vtk_paths = sorted(vtk_dir.glob("**/*.vtk"))

        console.print(vtk_paths)
        mesh_tc = cls.from_vtks(vtk_paths=vtk_paths, show=show)
        return mesh_tc

    @classmethod
    def from_vtks(cls, vtk_paths: Iterable[Path], show: bool = False) -> MeshTimecourse:
        """Reads mesh timecourse from multiple VTKs.

        :param vtk_paths: Path to VTK file.
        :param show: boolean flag to show information.
        """
        console.print(vtk_paths)
        mres_dict: Dict[str, MeshTimepoint] = {}

        # Nt: int = len(vtk_paths)
        # time = np.zeros(shape=(Nt,))

        for k, path in enumerate(vtk_paths):
            console.print(path)
            # time[k] = int(str(path).split(".")[-1][1:])

            mres: MeshTimepoint = cls.from_vtk(vtk_path=path, show=False)
            # mres_dict[path.name] = mres
            #
            # if k > 5:
            #     break

        # print(time)


def calculate_statistics(mp: MeshTimepoint, key: str):
    values: np.ndarray = mp.cell_data[key]
    if len(values.shape) != 1:
        raise ValueError("statistics can only be calculated on 1D data")

    # mean, sd, min, max (over time)


def calculate_zonated_statistics(mp: MeshTimepoint, key: str):
    pass


if __name__ == "__main__":
    from pm import EXAMPLE_VTK, DATA_DIR
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

