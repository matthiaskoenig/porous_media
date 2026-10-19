"""Helper functions for manipulation of geometries.

Reconstruction of geometries from partial simulations.
Helper function for composition and manipulation of mesh geometries.
"""

import os
from collections.abc import Mapping
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import meshio
import numpy as np
from numpy.typing import ArrayLike
from rich.progress import track

from porous_media.data.xdmf_tools import XDMFTimeSeriesWriter


def rotate_points(points: np.ndarray, angle: float) -> np.ndarray:
    """Rotate points clockwise by an angle in radians.

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
    # Exclude the origin and unwrap around a circular mean to handle +/- pi.
    xy_points = points[:, :2]
    radii = np.linalg.norm(xy_points, axis=1)
    xy_points = xy_points[radii > radii.max() * 1e-8]
    if not len(xy_points):
        raise ValueError("Sector points must extend beyond the origin")

    # Calculate angles with respect to the origin
    # angle of all points
    angles = np.arctan2(xy_points[:, 1], xy_points[:, 0])

    center = np.angle(np.mean(np.exp(1j * angles)))
    angles = center + np.angle(np.exp(1j * (angles - center)))
    result: float = float((angles.min() + angles.max()) / 2)
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


def sector_transforms(mirror_angle: float) -> tuple[np.ndarray, ...]:
    """Place six sectors counterclockwise, reflecting every other sector.

    Alternating reflections pair like radial boundaries. Seam points remain
    separate, so different patterns and resolutions retain their own values.
    Matrices act on column vectors; row arrays multiply their transpose.
    """
    reflection = np.array(
        [
            [np.cos(2 * mirror_angle), np.sin(2 * mirror_angle), 0.0],
            [np.sin(2 * mirror_angle), -np.cos(2 * mirror_angle), 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    transforms = []
    for sector in range(6):
        angle = sector * np.pi / 3
        rotation = np.array(
            [
                [np.cos(angle), -np.sin(angle), 0.0],
                [np.sin(angle), np.cos(angle), 0.0],
                [0.0, 0.0, 1.0],
            ]
        )
        transforms.append(rotation @ reflection if sector % 2 else rotation)
    return tuple(transforms)


def transform_cells(block: meshio.CellBlock, matrix: np.ndarray) -> np.ndarray:
    """Preserve the source cell orientation when reflecting a mesh."""
    if np.linalg.det(matrix) > 0:
        return np.asarray(block.data).copy()
    permutations = {
        "vertex": [0],
        "line": [1, 0],
        "triangle": [0, 2, 1],
        "quad": [0, 3, 2, 1],
        "tetra": [1, 0, 2, 3],
        "wedge": [0, 2, 1, 3, 5, 4],
        "hexahedron": [0, 3, 2, 1, 4, 7, 6, 5],
        "pyramid": [0, 3, 2, 1, 4],
    }
    if block.type not in permutations:
        raise ValueError(f"Unsupported reflected cell type: {block.type}")
    return np.asarray(block.data)[:, permutations[block.type]]


def transform_values(values: ArrayLike, matrix: np.ndarray) -> np.ndarray:
    """Transform Cartesian vectors and tensors, preserving scalar values.

    Full tensors use Q T Q.T. Unsupported component layouts are rejected rather
    than silently copying directional data into a different coordinate system.
    """
    array = np.asarray(values)
    dtype = np.result_type(array.dtype, np.float32)
    matrix = matrix.astype(dtype)
    if array.ndim == 1 or array.shape[1:] == (1,):
        return array.copy()
    if array.shape[1:] == (3,):
        return array @ matrix.T
    if array.shape[1:] in {(3, 3), (9,)}:
        tensor = array.reshape(-1, 3, 3)
        return (matrix @ tensor @ matrix.T).reshape(array.shape)
    raise ValueError(f"Unsupported field component shape: {array.shape[1:]}")


def transform_fields(
    data: Mapping[str, ArrayLike], matrix: np.ndarray
) -> dict[str, np.ndarray]:
    """Transform fields, including Cartesian coordinates stored as scalars."""
    transformed = {
        name: transform_values(values, matrix) for name, values in data.items()
    }
    names = ("position_x", "position_y", "position_z")
    if any(name in data for name in names):
        if not all(name in data for name in names):
            raise ValueError("Cartesian position fields must include x, y and z")
        positions = np.column_stack(
            [np.asarray(data[name]).reshape(-1) for name in names]
        )
        result = transform_values(positions, matrix)
        for index, name in enumerate(names):
            transformed[name] = result[:, index].reshape(np.asarray(data[name]).shape)
    return transformed


def create_lobulus_mesh(
    points: np.ndarray, cells: list[meshio.CellBlock]
) -> tuple[np.ndarray, list[meshio.CellBlock]]:
    """Reconstruct six alternating mirrored sectors without welding seams."""
    transforms = sector_transforms(find_reflection_angle_in_sixth(points))
    return np.concatenate([points @ matrix.T for matrix in transforms]), [
        meshio.CellBlock(
            block.type,
            np.concatenate(
                [
                    transform_cells(block, matrix) + len(points) * sector
                    for sector, matrix in enumerate(transforms)
                ]
            ),
        )
        for block in cells
    ]


def create_lobulus_data(
    point_data: dict[str, ArrayLike],
    cell_data: dict[str, list[ArrayLike]],
    transforms: tuple[np.ndarray, ...] | None = None,
) -> tuple[dict[str, ArrayLike], dict[str, list[ArrayLike]]]:
    """Replicate scalars and transform directional fields into each sector."""
    transforms = sector_transforms(0.0) if transforms is None else transforms
    point_sectors = [transform_fields(point_data, matrix) for matrix in transforms]
    new_point_data: dict[str, ArrayLike] = {
        name: np.concatenate([sector[name] for sector in point_sectors])
        for name in point_data
    }
    new_cell_data: dict[str, list[ArrayLike]] = {name: [] for name in cell_data}
    if cell_data:
        for block_index in range(len(next(iter(cell_data.values())))):
            block_data = {
                name: blocks[block_index] for name, blocks in cell_data.items()
            }
            sectors = [transform_fields(block_data, matrix) for matrix in transforms]
            for name in cell_data:
                new_cell_data[name].append(
                    np.concatenate([sector[name] for sector in sectors])
                )
    return new_point_data, new_cell_data


def reconstruct_lobulus_from_hexagon(xdmf_in: Path, xdmf_out: Path) -> None:
    """Construct lobulus data layer from sixth xdmf file."""
    if xdmf_in.resolve() == xdmf_out.resolve():
        raise ValueError("Input and output XDMF paths must differ")
    for path in (xdmf_out, xdmf_out.with_suffix(".h5")):
        if path.exists() or path.is_symlink():
            raise FileExistsError(f"Refusing to overwrite existing geometry: {path}")
    xdmf_out.parent.mkdir(parents=True, exist_ok=True)
    # Read the original sixth of the circle mesh
    with meshio.xdmf.TimeSeriesReader(xdmf_in) as reader:
        # create the mesh
        points, cells = reader.read_points_cells()
        new_points, new_cells = create_lobulus_mesh(points, cells)
        transforms = sector_transforms(find_reflection_angle_in_sixth(points))

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
                    transforms=transforms,
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
