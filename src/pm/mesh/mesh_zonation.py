"""Create zonated meshes for the analysis.
"""
from pathlib import Path
from typing import Callable, List, Dict

import meshio
import numpy as np

from pm.console import console
from pm.visualization.image_manipulation import merge_images
from pm.visualization.pyvista_visualization import visualize_lobulus_vtk


def create_zonated_mesh(mesh: meshio.Mesh, remove_point_data: bool=True, remove_cell_data: bool = True,
                        copy_mesh: bool = True) -> meshio.Mesh:
    """Calculates the distance from periportal and perivenous.

    Uses the cell_type variable for determining the position in the lobulus (periportal, perivenous).
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

    # calculate positions based on cell_type info
    cell_type: np.ndarray = m.cell_data["cell_type"][0]
    position = np.NaN * np.ones_like(cell_type)

    pp_cells: Dict[int, np.ndarray] = {}
    pv_cells: Dict[int, np.ndarray] = {}
    inner_cells: Dict[int, np.ndarray] = {}
    cell_block: meshio._mesh.CellBlock

    for (_, cell_block) in enumerate(m.cells):
        for kc, cell in enumerate(cell_block.data):
            # console.print(cell)

            # calculate the center of mass of cell
            count = 0
            center = np.zeros(shape=3)
            for kp in cell:
                center += m.points[kp]
                count += 1
            center = center / count

            # periportal cell
            if np.isclose(cell_type[kc], 1.0):
                position[kc] = 0
                pp_cells[kc] = center

            # pericentral cell
            elif np.isclose(cell_type[kc], 2.0):
                position[kc] = 1.0
                pv_cells[kc] = center

            # inner cell
            elif np.isclose(cell_type[kc], 0.0):
                position[kc] = np.NaN
                inner_cells[kc] = center

    # calculate pp-pv position for all inner cells
    for kc, center in inner_cells.items():
        # shortest distance periportal
        dpv = 1.0
        for center_pv in pv_cells.values():
            d = np.linalg.norm(center-center_pv)
            if d < dpv:
                dpv = d

        # shortest distance pericentral
        dpp = 1.0
        for center_pp in pp_cells.values():
            d = np.linalg.norm(center - center_pp)
            if d < dpp:
                dpp = d

        position[kc] = dpp/(dpv + dpp)

    console.print(f"{position=}")
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
        1: perivenous/pericentral
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

    # FIXME: use new patterns with constant amount

    @staticmethod
    def position(p: np.ndarray) -> np.ndarray:
        """Position information."""
        return p

    @staticmethod
    def constant(p: np.ndarray) -> np.ndarray:
        """Constant zonation."""
        return 0.5 * np.ones_like(p)

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
        data[p <= 0.2] = 1.0
        return data

    @staticmethod
    def sharp_pericentral(p: np.ndarray) -> np.ndarray:
        """Sharp pericentral pattern in [0, 1]."""
        data = np.zeros_like(p)
        data[p >= 0.8] = 1.0
        return data


def visualize_patterns(mesh: meshio.Mesh):
    # create vtk
    vtk_path = Path('../resources/zonation/mesh_zonation.vtk')
    mesh.write(vtk_path)

    output_path = Path("../../../results/zonation/raw_patterns/")
    output_path.mkdir(exist_ok=True)

    scalars = {}
    for key in mesh.cell_data:
        if key.startswith("pattern__") or key == "cell_type":
            pattern = key.split("__")[-1]
            data = mesh.cell_data[key][0]
            # new min, max
            dmin = data.min()
            dmax = data.max()
            # hardcoded for zonation patterns (FIXME: do on global data)
            # dmin = 0.0
            # dmax = 5.0

            scalars[key] = {
                "title": f"{pattern.upper()} [-]",
                # FIXME: better colormap
                "cmap": "RdBu",
                # "cmap": "Blues",
                "clim": (dmin, dmax),
            }

    # create raw images
    visualize_lobulus_vtk(
        vtk_path=vtk_path,
        scalars=scalars,
        output_dir=output_path
    )
    # combine images
    images: List[Path] = []
    for scalar in scalars:
        img_path = output_path / scalar / f"{vtk_path.stem}.png"
        images.append(img_path)

    image: Path = output_path / f"{vtk_path.stem}.png"
    merge_images(paths=images, direction="horizontal", output_path=image)
    console.print(f"Image created: {image}")


if __name__ == "__main__":
    # mesh
    vtk_path = Path(__file__).parent / "mesh_zonation.vtk"
    mesh: meshio.Mesh = meshio.read(vtk_path)
    m: meshio.Mesh = create_zonated_mesh(mesh)

    for f in [
        # ZonationPatterns.position,
        ZonationPatterns.constant,
        # ZonationPatterns.random,
        ZonationPatterns.linear_increase,
        ZonationPatterns.linear_decrease,
        # ZonationPatterns.exp_increase,
        ZonationPatterns.sharp_pericentral,
        # ZonationPatterns.exp_decrease,
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

    # mesh_tp: MeshTimepoint = MeshTimepoint.from_vtk(vtk_path=vtk_path, show=True)
    # console.rule(style="white")

    console.rule(title="Mesh Visualization", style="white")
    visualize_patterns(m2)
