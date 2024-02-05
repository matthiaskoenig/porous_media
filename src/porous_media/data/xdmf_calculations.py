"""Helper functions for calculating information on the Timecourse XDMF.

Includes things such as calculation of sum, ...
How to best handle this?
Provide some simple helpers for common operations, e.g.
- integral over the mesh (e.g. total necrosis, fat)
- relative fraction of the mesh (with max values or one)
- densities (normalization on volumes)

This should work on cell_data and point_data of the mesh

A lot of the information should be position dependent for analysis.
Support selections of positions; histogramm over position.

"""

from pathlib import Path
from typing import Dict, List, Optional, Tuple

import meshio
import numpy as np
import pandas as pd
import xarray as xr
from rich.progress import Progress, track

from porous_media import DATA_DIR, RESULTS_DIR
from porous_media.console import console
from porous_media.data.xdmf_tools import AttributeType, XDMFInfo


def mesh_datasets_from_xdmf(
    xdmf_path: Path,
    overwrite: bool = False,
) -> Tuple[xr.Dataset, xr.Dataset]:
    """Create the cell_data and point_data xarray DataSet.

    Here the actual mesh geometry is not important, but information is stored for
    cells and points. Files are serialized as netCDF4 for caching information.
    """
    xr_cells_path = xdmf_path.parent / f"{xdmf_path.stem}_cells.nc"
    xr_points_path = xdmf_path.parent / f"{xdmf_path.stem}_points.nc"
    if not overwrite and xr_cells_path.exists() and xr_points_path.exists():
        console.print(
            f"cells and points files exist: {xr_cells_path} | {xr_points_path}"
        )
        xr_cells: xr.Dataset = xr.open_dataset(xr_cells_path)
        xr_points: xr.Dataset = xr.open_dataset(xr_points_path)
        return xr_cells, xr_points

    xdmf_info: XDMFInfo = XDMFInfo.from_path(xdmf_path)
    # console.print(xdmf_info)

    dfs_cell_data: List[pd.DataFrame] = []
    dfs_point_data: List[pd.DataFrame] = []

    with meshio.xdmf.TimeSeriesReader(xdmf_path) as reader:
        data_cell: Dict[str, np.ndarray] = {}
        data_point: Dict[str, np.ndarray] = {}

        _, _ = reader.read_points_cells()
        tnum = reader.num_steps
        timepoints = np.ndarray(shape=(tnum,))
        for k in track(
            range(tnum),
            description=f"Create dataframes for mesh data '{xdmf_path}' ...",
        ):
            t, point_data, cell_data = reader.read_data(k)
            timepoints[k] = t

            # parse cell data
            for key, data in cell_data.items():
                # How to deal with vectors and tensors (FIXME)? (length of the vector)
                if xdmf_info.cell_data[key].attribute_type != AttributeType.SCALAR:
                    continue

                data = np.dstack(data)
                data = data.squeeze()
                data_cell[key] = data

            df_cell = pd.DataFrame(data_cell)
            dfs_cell_data.append(df_cell)

            # parse point data
            for key, data in point_data.items():
                if xdmf_info.point_data[key].attribute_type != AttributeType.SCALAR:
                    continue

                data = np.dstack(data)
                data = data.squeeze()
                data_point[key] = data

            # FIXME: directly create the xarray datasets;
            df_point = pd.DataFrame(data_point)
            dfs_point_data.append(df_point)

    with Progress() as progress:
        _ = progress.add_task("Create xarray datasets ...", total=None)

        # (cell, time) xarray Dataset
        xr_cells = xr.concat(
            [df.to_xarray() for df in dfs_cell_data],
            dim="time",
        )
        xr_cells = xr_cells.rename_dims(dims_dict={"index": "cell"})
        xr_cells = xr_cells.assign_coords(coords={"time": timepoints})
        # serialize to netCDF4
        xr_cells.to_netcdf(xr_cells_path)

        # (point, time) xarray Dataset
        xr_points = xr.concat(
            [df.to_xarray() for df in dfs_point_data],
            dim="time",
        )
        xr_points = xr_points.rename_dims(dims_dict={"index": "point"})
        xr_points = xr_points.assign_coords(coords={"time": timepoints})
        # serialize to netCDF4
        xr_points.to_netcdf(xr_points_path)

    return xr_cells, xr_points


if __name__ == "__main__":
    console.rule(title="XDMF calculations", style="white")

    # prepare data for analysis
    xdmf_path = (
        Path("/home/mkoenig/git/porous_media/data/spt/2023-12-05") / "sim001.xdmf"
    )
    xr_cells, xr_points = mesh_datasets_from_xdmf(xdmf_path, overwrite=False)

    console.rule(align="left", title="cell_data", style="white")
    console.print(xr_cells)
    console.rule(align="left", title="point_data", style="white")
    console.print(xr_points)
