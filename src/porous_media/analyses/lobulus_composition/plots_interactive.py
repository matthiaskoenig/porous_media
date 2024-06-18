# Simple test that all variables are read
import os
from pathlib import Path
from typing import List, Tuple
import shutil
import meshio
import numpy as np

from porous_media import DATA_DIR
from porous_media.console import console
from porous_media.data.xdmf_tools import XDMFInfo, AttributeType
from porous_media.visualization.pyvista_visualization import DataLayer, xdmf_to_mesh, \
    visualize_interactive, VisualizationSettings
from rich.progress import track


def rotate_points(points: np.ndarray, angle: float) -> np.ndarray:
    """Rotate points clockwise by a certain angle in degrees.

    In linear algebra, a rotation matrix is a transformation matrix that is used
    to perform a rotation in Euclidean space.

    angle in radians between 0 and 2*pi
    """
    rotation_matrix = np.array([
        [np.cos(angle), -np.sin(angle), 0.0],
        [np.sin(angle), np.cos(angle), 0.0],
        [0.0, 0.0, 1.0],
    ])

    return np.dot(points, rotation_matrix)


def create_lobulus_mesh(points: np.ndarray, cells: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Create the lobulus from the sixth element."""
    # Duplicate and rotate the mesh five times
    angles = 2 * np.pi / 6 * np.linspace(start=1, stop=5, num=5)

    all_points = [points]
    all_cells = [cells]

    # create points
    n_points = len(points)
    for k, angle in enumerate(angles):
        # Rotate only x and y
        rotated_points = rotate_points(points, angle)
        all_points.append(rotated_points)

        # Adjust cell connectivity
        # find points and connect
        offset = n_points * (k + 1)
        new_cells = []

        for cell_block in cells:
            new_data = cell_block.data + offset
            new_cells.append(meshio.CellBlock(cell_block.type, new_data))

        all_cells.extend([new_cells])

    full_circle_points = np.vstack(all_points)

    # create cells
    full_circle_cells = []
    cell_types = set(cell_block.type for cell_block in cells)

    for cell_type in cell_types:
        cells_for_type = []
        cell_block: meshio.CellBlock
        for cell_block_list in all_cells:
            for cell_block in cell_block_list:
                if cell_block.type == cell_type:
                    cells_for_type.append(cell_block.data)

        full_circle_cells.append(
            meshio.CellBlock(cell_type, np.vstack(cells_for_type))
        )

    return full_circle_points, full_circle_cells


def create_lobulus_data(point_data: np.ndarray, cell_data: List[np.ndarray]) -> Tuple:
    """Create the lobulus data from the sixth element."""
    new_point_data = {}
    for key, data in point_data.items():
        new_point_data[key] = np.vstack([data for _ in range(6)])

    new_cell_data = {}
    for key, data_list in cell_data.items():
        new_data_list = []
        for data in data_list:
            new_data_list.append(np.vstack([data for _ in range(6)]))

        new_cell_data[key] = new_data_list

    return new_point_data, new_cell_data


def reconstruct_lobulus_from_hexagon(xdmf_in: Path, xdmf_out: Path) -> DataLayer:
    """Construct lobulus data layer from hexagonal xdmf file."""

    # Read the original sixth of the circle mesh
    with meshio.xdmf.TimeSeriesReader(xdmf_in) as reader:
        points, cells = reader.read_points_cells()
        new_points, new_cells = create_lobulus_mesh(points, cells)

        with meshio.xdmf.TimeSeriesWriter(xdmf_out, data_format="HDF") as writer:
            writer.write_points_cells(new_points, new_cells)

            # process all the data
            for k in track(range(reader.num_steps), description="Process data"):
                t, point_data, cell_data = reader.read_data(k)
                new_point_data, new_cell_data = create_lobulus_data(
                    point_data=point_data,
                    cell_data=cell_data,
                )
                writer.write_data(t, point_data=new_point_data, cell_data=new_cell_data)


    # Fix incorrect *.h5 path
    # https://github.com/nschloe/meshio/pull/1358
    h5_path = xdmf_out.parent / f"{xdmf_out.stem}.h5"
    if h5_path.exists():
        os.remove(h5_path)
    shutil.move(f"{xdmf_out.stem}.h5", str(xdmf_out.parent))


xdmf_path = DATA_DIR / "lobulus_composition" / "sim003.xdmf"
xdmf_lobulus_path = DATA_DIR / "lobulus_composition" / "sim003_lobulus.xdmf"
reconstruct_lobulus_from_hexagon(xdmf_in=xdmf_path, xdmf_out=xdmf_lobulus_path)

xdmf_info: XDMFInfo = XDMFInfo.from_path(xdmf_lobulus_path)
console.print(xdmf_info)


# --- Visualization

data_layers: List[DataLayer] = [
    DataLayer(
        sid="displacement",
        title="displacement",
        colormap="magma",
        viz_type=AttributeType.VECTOR,
    ),
    DataLayer(
        sid="fluid_flux_TPM",
        title="fluid flux [m/s]",
        colormap="magma",
        viz_type=AttributeType.VECTOR,
    ),
    DataLayer(
        sid="pressure",
        title="Pressure",
        colormap="magma",
        viz_type=AttributeType.SCALAR,
    ),
    DataLayer(
        sid="effective_fluid_pressure_TPM",
        title="Effective fluid pressure [Pa]",
        colormap="magma",
        viz_type=AttributeType.SCALAR,
    ),
    DataLayer(
        sid="rr_(S)",
        title="Substrate S [mM]",
        colormap="magma",
        viz_type=AttributeType.SCALAR,
    ),
]
data_layers_dict = {dl.sid: dl for dl in data_layers}

mesh = xdmf_to_mesh(xdmf_path, k=1)
mesh = xdmf_to_mesh(xdmf_lobulus_path, k=1)
console.print(mesh)
visualize_interactive(
    mesh, data_layer=data_layers_dict["rr_(S)"], visualization_settings=VisualizationSettings()
)

