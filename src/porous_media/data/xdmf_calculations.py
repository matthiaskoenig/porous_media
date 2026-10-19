"""Convert scalar XDMF time series into in-memory xarray datasets."""

import os
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from multiprocessing import get_context
from pathlib import Path

import meshio
import numpy as np
import xarray as xr


def _scalar_fields(data: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    """Select scalar arrays without squeezing away single-element dimensions."""
    return {
        name: values.reshape(-1)
        for name, values in data.items()
        if values.ndim == 1 or (values.ndim == 2 and values.shape[1] == 1)
    }


def _store_step(
    target: dict[str, np.ndarray],
    fields: dict[str, np.ndarray],
    step: int,
    num_steps: int,
    num_items: int,
) -> None:
    """Fill preallocated arrays while checking schema and preserving promoted types."""
    if step and target.keys() != fields.keys():
        raise ValueError("Scalar fields must remain consistent across time steps")
    for name, values in fields.items():
        if values.size != num_items:
            raise ValueError(f"Scalar field {name!r} has an inconsistent size")
        if step == 0:
            target[name] = np.empty((num_steps, num_items), dtype=values.dtype)
        elif not np.can_cast(values.dtype, target[name].dtype, casting="safe"):
            target[name] = target[name].astype(
                np.result_type(target[name].dtype, values.dtype)
            )
        target[name][step] = values


def _cache_mesh_datasets(xdmf_path: Path, overwrite: bool) -> None:
    """Write caches without sending large datasets back through a process pool."""
    mesh_datasets_from_xdmf(xdmf_path, overwrite=overwrite)


def mesh_datasets_from_xdmf_dir(
    xdmf_dir: Path, overwrite: bool = True, *, max_workers: int | None = None
) -> None:
    """Cache scalar datasets using a bounded number of worker processes.

    Args:
        xdmf_dir: Directory containing XDMF time series.
        overwrite: Whether to regenerate existing caches.
        max_workers: Positive worker limit. Defaults to available CPUs minus one,
            bounded by the number of files. One worker runs without a process pool.

    Raises:
        ValueError: The worker limit is not positive.
    """
    if max_workers is not None and max_workers < 1:
        raise ValueError("max_workers must be positive")
    paths = sorted(xdmf_dir.glob("*.xdmf"))
    if not paths:
        return
    workers = min(len(paths), max_workers or max(1, (os.process_cpu_count() or 1) - 1))
    if workers == 1:
        for path in paths:
            _cache_mesh_datasets(path, overwrite)
        return
    with ProcessPoolExecutor(
        max_workers=workers, mp_context=get_context("spawn")
    ) as executor:
        list(executor.map(partial(_cache_mesh_datasets, overwrite=overwrite), paths))


def mesh_datasets_from_xdmf(
    xdmf_path: Path, overwrite: bool = False
) -> tuple[xr.Dataset, xr.Dataset]:
    """Return scalar cell and point datasets and cache them as NetCDF files.

    Scalars from all cell blocks are concatenated in mesh order. Vector and tensor
    fields are omitted. Datasets use ``(time, cell)`` or ``(time, point)`` dimensions,
    including single-item and empty meshes. Time steps must share a scalar schema.
    Both newly generated and cached datasets are fully loaded; no open file handles
    are returned. Use ``overwrite=True`` after modifying source data.

    Args:
        xdmf_path: Fixed-mesh XDMF time series.
        overwrite: Regenerate the adjacent ``_cells.nc`` and ``_points.nc`` caches.

    Returns:
        Cell and point datasets, in that order.

    Raises:
        ValueError: The time series is empty or its scalar schema changes.
    """
    cells_path = xdmf_path.with_name(f"{xdmf_path.stem}_cells.nc")
    points_path = xdmf_path.with_name(f"{xdmf_path.stem}_points.nc")
    if not overwrite and cells_path.is_file() and points_path.is_file():
        return xr.load_dataset(cells_path), xr.load_dataset(points_path)

    cells_data: dict[str, np.ndarray] = {}
    points_data: dict[str, np.ndarray] = {}
    with meshio.xdmf.TimeSeriesReader(xdmf_path) as reader:
        points, cells = reader.read_points_cells()
        if reader.num_steps == 0:
            raise ValueError("XDMF time series must contain at least one time step")
        times = np.empty(reader.num_steps)
        num_cells = sum(len(block) for block in cells)
        for step in range(reader.num_steps):
            time, point_data, cell_data = reader.read_data(step)
            times[step] = time
            cell_fields = {
                name: np.concatenate(blocks, axis=0)
                for name, blocks in cell_data.items()
                if blocks
            }
            _store_step(
                cells_data,
                _scalar_fields(cell_fields),
                step,
                reader.num_steps,
                num_cells,
            )
            _store_step(
                points_data,
                _scalar_fields(point_data),
                step,
                reader.num_steps,
                len(points),
            )

    xr_cells = xr.Dataset(
        {name: (("time", "cell"), values) for name, values in cells_data.items()},
        coords={"time": times, "cell": np.arange(num_cells)},
    )
    xr_points = xr.Dataset(
        {name: (("time", "point"), values) for name, values in points_data.items()},
        coords={"time": times, "point": np.arange(len(points))},
    )
    xr_cells.to_netcdf(cells_path, engine="scipy")
    xr_points.to_netcdf(points_path, engine="scipy")
    return xr_cells, xr_points
