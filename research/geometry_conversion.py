"""Helper functions for manipulation of geometries.

Reconstruction of geometries from partial simulations.
Helper function for composition and manipulation of mesh geometries.
"""

import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import meshio
import numpy as np
from numpy.typing import ArrayLike
from rich.progress import track

from porous_media.data.xdmf_tools import XDMFTimeSeriesWriter


def rotate_points(points: np.ndarray, angle: float) -> np.ndarray:
    """Rotate points clockwise by a certain angle in degrees.

    In linear algebra, a rotation matrix is a transformation matrix that is used
    to perform a rotation in Euclidean space.

    angle in radians between 0 and 2*pi
    """
    rotation_matrix = np.array(
        [
            [np.cos(angle), -np.sin(angle), 0.0],
            [np.sin(angle), np.cos(angle), 0.0],
            [0.0, 0.0, 1.0],
        ]
    )

    result: np.ndarray = np.dot(points, rotation_matrix)
    return result


def find_reflection_angle_in_sixth(points: np.ndarray) -> float:
    """Find the reflection angle from points through the origin."""
    # Extract x and y coordinates
    xy_points = points[:, :2]

    # Calculate angles with respect to the origin
    # angle of all points
    angles = np.arctan2(xy_points[:, 1], xy_points[:, 0])

    # find mean angle for flipping
    result: float = (np.nanmin(angles) + np.nanmax(angles)) / 2
    return result


def reflect_points(points: np.ndarray, angle: float) -> np.ndarray:
    """Reflect points at mirror line defined by angle from origin."""
    # Construct the 3D reflection matrix for the given angle
    reflection_matrix = np.array(
        [
            [np.cos(2 * angle), np.sin(2 * angle), 0],
            [np.sin(2 * angle), -np.cos(2 * angle), 0],
            [0, 0, 1],
        ]
    )

    # Apply the reflection matrix to each point
    result: np.ndarray = np.dot(points, reflection_matrix.T)
    return result


def minimal_distance(points: np.ndarray) -> float:
    """Calculate the minimal distance between points."""
    # Calculate the pairwise distance matrix
    distances = np.linalg.norm(points[:, np.newaxis] - points, axis=2)

    # Set the diagonal to infinity to ignore zero distances between the same points
    np.fill_diagonal(distances, np.inf)

    # Find the minimum distance
    min_distance: float = np.min(distances)

    return min_distance


def unique_with_tolerance(
    points: np.ndarray, tol: float
) -> tuple[np.ndarray, np.ndarray]:
    """Find the unique points within a given tolerance."""
    if not np.isfinite(tol) or tol <= 0:
        raise ValueError("tol must be a positive finite number")

    unique_points = []
    unique_inverse = np.full(points.shape[0], -1, dtype=int)
    current_label = 0
    for i, point in enumerate(points):
        if unique_inverse[i] == -1:
            distances = np.linalg.norm(points - point, axis=1)
            close_points = (distances < tol) & (unique_inverse == -1)
            unique_inverse[close_points] = current_label
            unique_points.append(points[close_points].mean(axis=0))
            current_label += 1

    unique_points = np.array(unique_points) if unique_points else points[:0].copy()
    return unique_points, unique_inverse


def remove_duplicate_points(
    mesh: meshio.Mesh, unique_points: np.ndarray, unique_indices: np.ndarray
) -> meshio.Mesh:
    """Remap cells and keep data from the first point in each duplicate group."""
    # Step 2: Update cell data to use new point indices
    new_cells: list[tuple[str, ArrayLike] | meshio.CellBlock] = []
    for cell_block in mesh.cells:
        new_cell_data = unique_indices[cell_block.data]
        new_cells.append(meshio.CellBlock(cell_block.type, new_cell_data))
    # Step 3: Filter point data to keep only the rows corresponding to unique points
    _, first_indices = np.unique(unique_indices, return_index=True)
    new_point_data: dict[str, ArrayLike] = {}
    for key, data in mesh.point_data.items():
        new_point_data[key] = np.asarray(data)[first_indices]
    # Step 4: Create a new mesh with unique points and updated cells
    return meshio.Mesh(
        points=unique_points,
        cells=new_cells,
        point_data=new_point_data,
        cell_data=mesh.cell_data,
        field_data=mesh.field_data,
    )


def create_lobulus_mesh(
    points: np.ndarray, cells: list[meshio.CellBlock]
) -> tuple[np.ndarray, list[meshio.CellBlock]]:
    """Create the lobulus from the sixth element.

    Requires to copy all the data associated with the geometry.
    """
    # Duplicate and rotate the mesh five times
    angles = 2 * np.pi / 6 * np.linspace(start=1, stop=5, num=5)

    all_points = [points]

    # create points
    n_points = len(points)
    mirror_angle = find_reflection_angle_in_sixth(points)
    for angle in angles:
        # mirror for matching points
        reflected_points = reflect_points(points, mirror_angle)

        # Rotate only x and y
        rotated_points = rotate_points(reflected_points, angle)
        all_points.append(rotated_points)

    full_circle_points = np.vstack(all_points)

    # Keep the original block order so connectivity and cell data stay aligned.
    full_circle_cells = [
        meshio.CellBlock(
            block.type,
            np.vstack([np.asarray(block.data) + n_points * k for k in range(6)]),
        )
        for block in cells
    ]

    # FIXME
    # # remove duplicate points
    # mesh = meshio.Mesh(
    #     points=points, cells=cells, point_data=point_data, cell_data=cell_data
    # )
    # # Remove duplicate points and update cells
    # modified_mesh = remove_duplicate_points(mesh, unique_points, unique_indices)

    return full_circle_points, full_circle_cells


def create_lobulus_data(
    point_data: dict[str, ArrayLike], cell_data: dict[str, list[ArrayLike]]
) -> tuple[dict[str, ArrayLike], dict[str, list[ArrayLike]]]:
    """Create the lobulus data from the sixth element."""
    new_point_data: dict[str, ArrayLike] = {}
    for key, data in point_data.items():
        new_point_data[key] = np.concatenate([np.asarray(data)] * 6)

    new_cell_data: dict[str, list[ArrayLike]] = {}
    for key, data_list in cell_data.items():
        new_cell_data[key] = [
            np.concatenate([np.asarray(data)] * 6) for data in data_list
        ]

    return new_point_data, new_cell_data


def reconstruct_lobulus_from_hexagon(xdmf_in: Path, xdmf_out: Path) -> None:
    """Construct lobulus data layer from sixth xdmf file."""
    xdmf_out.parent.mkdir(parents=True, exist_ok=True)
    # Read the original sixth of the circle mesh
    with meshio.xdmf.TimeSeriesReader(xdmf_in) as reader:
        # create the mesh
        points, cells = reader.read_points_cells()
        new_points, new_cells = create_lobulus_mesh(points, cells)

        with XDMFTimeSeriesWriter(xdmf_out, data_format="HDF") as writer:
            writer.write_points_cells(
                new_points, [(block.type, block.data) for block in new_cells]
            )

            # process all the data
            for k in track(range(reader.num_steps), description="Process data"):
                t, point_data, cell_data = reader.read_data(k)
                new_point_data, new_cell_data = create_lobulus_data(
                    point_data=point_data,
                    cell_data=cell_data,
                )
                writer.write_data(t, point_data=new_point_data, cell_data=new_cell_data)


def reconstruct_lobulus_from_hexagon_dir(xdmf_in_dir: Path, xdmf_out_dir: Path) -> None:
    """Parallel reconstruction of lobulus."""
    xdmf_out_dir.mkdir(parents=True, exist_ok=True)

    files = [
        (xdmf_in, Path(xdmf_out_dir / xdmf_in.name))
        for xdmf_in in xdmf_in_dir.glob("*.xdmf")
    ]
    # parallel processing
    with ProcessPoolExecutor(max_workers=os.cpu_count()) as executor:
        futures = {
            executor.submit(
                reconstruct_lobulus_from_hexagon, xdmf_in, xdmf_out
            ): xdmf_in
            for xdmf_in, xdmf_out in files
        }
        for future in as_completed(futures):
            _ = future.result()

    # serial
    # for xdmf_in, xdmf_out in files:
    #     reconstruct_lobulus_from_hexagon(xdmf_in=xdmf_in, xdmf_out=xdmf_out)
