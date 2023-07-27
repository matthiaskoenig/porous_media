"""Create zonated meshes for the analysis.
"""
from pathlib import Path
from typing import Callable, List

import meshio
import numpy as np

from pm.console import console
from pm.mesh.mesh_tools import MeshTimepoint
from pm.visualization.image_manipulation import merge_images
from pm.visualization.pyvista_visualization import visualize_lobulus_vtk


def create_zonated_mesh(mesh: meshio.Mesh, remove_point_data: bool=True, remove_cell_data: bool = True, copy_mesh: bool = True) -> meshio.Mesh:
    """Calculates the distance from periportal and perivenous.

    Uses the cell_type variable for determining the distances.
    'cell_type':
        0: internal node
        1: periportal (inflow)
        2: perivenous (outflow)
    """

    # clone mesh
    if copy_mesh:
        m: meshio.Mesh = mesh.copy()
    else:
        m = mesh

    # check for variables
    if "cell_type" not in m.cell_data:
        raise IOError("'cell_type' required in cell_data for zonation patterns")

    # remove cell data
    if remove_cell_data:
        m.cell_data = {
            "cell_type": m.cell_data["cell_type"]
        }

    # remove point data
    if remove_point_data:
        m.point_data = {}

    # TODO: calculate the mesh point distance from periportal and perivenous
    console.print(f"Cells {m.cells}")
    for (k, cell_block) in enumerate(m.cells):
        console.print(cell_block)
        console.print(type(cell_block))

    # periportal cells
    # perivenous cells

    # distance of cells
    # 626 points

    # store the position
    cell_type: np.ndarray = m.cell_data["cell_type"][0]
    position = np.random.rand(*cell_type.shape)
    console.print(f"{position}")
    m.cell_data["position"] = [position]
    return m


def add_zonated_variable(
    mesh: meshio.Mesh,
    variable_id: str,
    f_zonation: Callable
) -> meshio.Mesh:
    """Uses the position variable in [0,1] to add a zonation variable to the mesh.

    :param variable_id: identifier in cell_data
    :param f_zonation: function to calulcate zonation

    position:
        0: periportal
        1: perivenous
    """
    # check for variables
    m = mesh
    if "position" not in m.cell_data:
        raise IOError("'position' required in calculate zonation patterns, 'create_zonated_mesh' first.")

    position: np.ndarray = m.cell_data["position"][0]
    data = f_zonation(position)
    m.cell_data[variable_id] = [data]
    return m


class ZonationPatterns:
    """Definition of standard zonation patterns."""

    @staticmethod
    def constant(p: np.ndarray) -> np.ndarray:
        """Constant zonation."""
        return np.ones_like(p)

    @staticmethod
    def random(p: np.ndarray) -> np.ndarray:
        """Random zonation in [0, 1]."""
        return np.random.rand(*p.shape)

    @staticmethod
    def linear_increase(p: np.ndarray) -> np.ndarray:
        """Linear increasing pattern in [0, 1]."""
        return p

    @staticmethod
    def linear_decrease(p: np.ndarray) -> np.ndarray:
        """Linear decreasing pattern in [0, 1]."""
        return 1.0 - ZonationPatterns.linear_increase(p)

    @staticmethod
    def exp_increase(p: np.ndarray) -> np.ndarray:
        """Exponential increasing pattern in [0, 1]."""
        return (np.exp(p) - 1.0)/(np.exp(1.0) - 1.0)

    @staticmethod
    def exp_decrease(p: np.ndarray) -> np.ndarray:
        """Exponential decreasing pattern in [0, 1]."""
        return 1.0 - ZonationPatterns.exp_increase(p)

    @staticmethod
    def sharp_periportal(p: np.ndarray) -> np.ndarray:
        """Sharp periportal pattern in [0, 1]."""
        data = np.zeros_like(p)
        data[p<=0.2] = 1.0
        return data

    @staticmethod
    def sharp_pericentral(p: np.ndarray) -> np.ndarray:
        """Sharp pericentral pattern in [0, 1]."""
        data = np.zeros_like(p)
        data[p>=0.8] = 1.0
        return data

def visualize_patterns(mesh: meshio.Mesh):
    # create vtk
    vtk_path = Path('mesh_zonation.vtk')
    mesh.write(vtk_path)

    output_path = Path("./raw_patterns/")
    output_path.mkdir(exist_ok=True)

    scalars = {}
    for key in mesh.cell_data:
        if key.startswith("pattern__"):
            pattern = key.split("__")[-1]
            data = mesh.cell_data[key][0]
            # new min, max
            dmin = data.min()
            dmax = data.max()
            scalars[key] = {
                "title": f"Zonation {pattern.upper()} [-]",
                "cmap": "RdBu",
                "clim": (dmin, dmax),
            }

    # create raw images
    visualize_lobulus_vtk(
        vtk_path=vtk_path,
        scalars=scalars,
        output_dir=output_path
    )
    # combine images
    row: List[Path] = []
    for scalar in scalars:
        img_path = output_path / scalar / f"{vtk_path.stem}.png"
        row.append(img_path)

    row_image: Path = output_path / f"{vtk_path.stem}_row.png"
    merge_images(paths=row, direction="horizontal", output_path=row_image)


if __name__ == "__main__":
    vtk_path = Path(__file__).parent / "mesh_zonation.vtk"
    mesh: meshio.Mesh = meshio.read(vtk_path)
    m: meshio.Mesh = create_zonated_mesh(mesh)

    for f in [
        ZonationPatterns.constant,
        ZonationPatterns.random,
        ZonationPatterns.linear_increase,
        ZonationPatterns.linear_decrease,
        ZonationPatterns.exp_increase,
        ZonationPatterns.exp_decrease,
        ZonationPatterns.sharp_pericentral,
        ZonationPatterns.sharp_periportal,
    ]:
        add_zonated_variable(
            mesh=m,
            variable_id=f"pattern__{f.__name__}",
            f_zonation=f,
        )

    console.rule(title="Mesh Serialization", style="white")
    console.print(f"{m=}")
    m.write(Path(__file__).parent / "mesh_zonation.xdmf")
    m2 = meshio.read(Path(__file__).parent / "mesh_zonation.xdmf")
    console.print(f"{m2=}")

    console.rule(style="white")
    # mesh_tp: MeshTimepoint = MeshTimepoint.from_vtk(vtk_path=vtk_path, show=True)
    # console.rule(style="white")

    console.rule(title="Mesh Visualization", style="white")
    visualize_patterns(m2)
