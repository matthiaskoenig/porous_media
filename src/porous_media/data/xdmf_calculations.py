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
from matplotlib import pyplot as plt
from pathlib import Path
from typing import Tuple, Dict, List

import meshio
import numpy as np
import pandas as pd
import xarray as xr
from rich.progress import track

from porous_media import DATA_DIR, RESULTS_DIR
from porous_media.console import console
from porous_media.data.xdmf_tools import XDMFInfo, AttributeType


def create_mesh_dataframes(xdmf_path: Path) -> Tuple[xr.Dataset, xr.Dataset]:
    """Create the cell_data and point_data DataFrame.

    Here the actual mesh geometry does not matter any more, but we integrate over all the elements.

    depending on the timepoints this will be a pretty large dataframe.

    """
    xdmf_info: XDMFInfo = XDMFInfo.from_path(xdmf_path)
    console.print(xdmf_info)

    dfs_cell_data: List[pd.DataFrame] = []
    dfs_point_data: List[pd.DataFrame] = []

    with meshio.xdmf.TimeSeriesReader(xdmf_path) as reader:

        data_cell: Dict[str, np.ndarray] = {}
        data_point: Dict[str, np.ndarray] = {}

        _, _ = reader.read_points_cells()
        tnum = reader.num_steps
        timepoints = np.ndarray(shape=(tnum,))
        for k in track(range(tnum), description=f"Creating DataFrames ..."):
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

    # generate xarrays

    # (cell, time) xarray Dataset
    xr_cells: xr.Dataset = xr.concat(
        [df.to_xarray() for df in dfs_cell_data],
        dim="time",
    )
    xr_cells = xr_cells.rename_dims(dims_dict={"index": "cell"})
    xr_cells = xr_cells.assign_coords(coords={"time": timepoints})

    # (point, time) xarray Dataset
    xr_points: xr.Dataset = xr.concat(
        [df.to_xarray() for df in dfs_point_data],
        dim="time",
    )
    xr_points = xr_points.rename_dims(dims_dict={"index": "point"})
    xr_points = xr_points.assign_coords(coords={"time": timepoints})

    return xr_cells, xr_points


def plot_necrosis_over_time(xr_cells_dict: Dict[str, xr.Dataset]):
    """Plot necrosis over time."""
    console.rule(title="necrosis calculation", style="white")
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(10, 10))

    # [1] necrosis fraction ~ time
    axes[0, 0].set_xlabel("time [hr]")
    axes[0, 0].set_ylabel("necrosis fraction [%]")

    for label, xr_cells in xr_cells_dict.items():
        # calculate necrosis fraction (sum/count)
        # FIXME: calculate and add the cell volumes for proper normalization
        necrosis = xr_cells.rr_necrosis
        necrosis_fraction = necrosis.sum(dim="cell")/necrosis.count(dim="cell")
        console.print(f"{necrosis_fraction=}")

        axes[0, 0].plot(
            # convert to hr and percent
            necrosis_fraction.time/60/60, necrosis_fraction * 100,
            label=label,
            linestyle="-",
            marker="o",
            markeredgecolor="black",
        )
    axes[0, 0].legend()

    # necrosis ~ position

    plt.show()


if __name__ == "__main__":
    console.rule(title="XDMF calculations", style="white")
    # interpolated dataframe for zonation patterns
    xdmf_paths = [RESULTS_DIR / "spt_zonation_patterns_new" / "10_28800.0" / f"simulation_pattern{k}_interpolated.xdmf" for k in range(5)]
    labels = [
        "Constant",
        "Linear increase",
        "Linear decrease",
        "Sharp pericentral",
        "Sharp periportal",
    ]

    # example dataset
    xdmf_path = xdmf_paths[1]
    xr_cells, xr_points = create_mesh_dataframes(xdmf_path)
    console.rule(align="left", title="cell_data", style="white")
    console.print(xr_cells)
    console.rule(align="left", title="point_data", style="white")
    console.print(xr_points)

    # calculate all xarray Datasets
    xr_cells_dict: Dict[str, xr.Dataset] = {}
    for k, xdmf_path in enumerate(xdmf_paths):
        label = labels[k]
        xr_cells, xr_points = create_mesh_dataframes(xdmf_path)
        xr_cells_dict[label] = xr_cells

    # calculate the necrosis area for all the simulations

    plot_necrosis_over_time(
        xr_cells_dict=xr_cells_dict
    )

