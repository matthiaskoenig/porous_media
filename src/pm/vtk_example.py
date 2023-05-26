from pathlib import Path
from typing import Dict, Tuple

import numpy as np

from pm.console import console
import meshio
import pyvista as pv

from dataclasses import dataclass, field


@dataclass
class MeshResults:
    mesh: meshio.Mesh
    cell_data: Dict[str, np.ndarray] = field(repr=False)
    cell_data_shape: Dict[str, Tuple]
    point_data: Dict[str, np.ndarray] = field(repr=False)
    point_data_shape: Dict[str, Tuple]


def read_mesh_data_from_vtk(vtk_path: Path) -> MeshResults:
    """Reads mesh information from VTK."""
    console.rule(align="left", style="white")
    console.print(vtk_path)
    mesh: meshio.Mesh = meshio.read(vtk_path)
    # data can be accessed from the mesh via
    # print("points:", mesh.points)
    # print("cells:", mesh.cells)
    # print("data:", mesh.point_data.keys())
    # print("data:", mesh.cell_data.keys())

    cell_data: Dict[str, np.ndarray] = {}
    cell_data_shape: Dict[str, np.ndarray] = {}
    for key in mesh.cell_data.keys():
        data = process_cell_data(key=key, mesh=mesh)
        cell_data[key] = data
        cell_data_shape[key] = data.shape

    point_data: Dict[str, np.ndarray] = {}
    point_data_shape: Dict[str, np.ndarray] = {}
    for key in mesh.point_data.keys():
        data = process_point_data(key=key, mesh=mesh)
        point_data[key] = data
        point_data_shape[key] = data.shape

    # console.print(f"{key: <20}: {type(data)}, {data.shape}")

    mesh_results = MeshResults(
        mesh=mesh,
        cell_data=cell_data,
        cell_data_shape=cell_data_shape,
        point_data=point_data,
        point_data_shape=point_data_shape,
    )
    console.print(mesh_results)
    return mesh_results


def process_cell_data(key: str, mesh: meshio.Mesh) -> np.ndarray:
    """Process data from the mesh."""
    data = mesh.cell_data[key]
    data = np.dstack(data)
    return data.squeeze()


def process_point_data(key: str, mesh: meshio.Mesh) -> np.ndarray:
    """Process data from the mesh."""
    data = mesh.point_data[key]
    data = np.dstack(data)
    return data.squeeze()


def visualize_mesh_pyvista(mesh: meshio.Mesh) -> None:
    """Blocking visualization with pyvista."""
    # simply pass the numpy points to the PolyData constructor
    # https://docs.pyvista.org/version/stable/examples/00-load/create-tri-surface.html#sphx-glr-examples-00-load-create-tri-surface-py
    cloud = pv.PolyData(mesh.points)
    cloud.plot(point_size=15)


if __name__ == "__main__":
    from pm import EXAMPLE_VTK
    mesh_results: MeshResults = read_mesh_data_from_vtk(EXAMPLE_VTK)

    # visualize_mesh_pyvista(mesh_results.mesh)




