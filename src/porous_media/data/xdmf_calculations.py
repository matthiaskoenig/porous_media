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
from typing import Tuple, Dict

import meshio
import numpy as np
import pandas as pd
from rich.progress import track

from porous_media import DATA_DIR
from porous_media.console import console
from porous_media.data.xdmf_tools import XDMFInformation


def create_mesh_dataframes(xdmf_path: Path) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Create the cell_data and point_data DataFrame.

    Here the actual mesh geometry does not matter any more, but we integrate over all the elements.

    depending on the timepoints this will be a pretty large dataframe.

    """
    df_cells = None
    df_points = None

    data_cells: Dict[str, np.ndarray] = {}
    data_points: Dict[str, np.ndarray] = {}

    with meshio.xdmf.TimeSeriesReader(xdmf_path) as reader:
        points, cells = reader.read_points_cells()

        tnum = reader.num_steps
        for k in track(range(tnum), description=f"Creating DataFrames ..."):
            t, point_data, cell_data = reader.read_data(k)

            # parse cell data
            for key, data in cell_data.items():
                console.print(key)
                data = np.dstack(data)
                data = data.squeeze()
                console.print(data.shape)
                cell_data[key] = data


            # parse point data
            for key, data in point_data.items():
                console.print(key)
                data = np.dstack(data)
                data = data.squeeze()
                console.print(data.shape)
                point_data[key] = data

            # FIXME: handle the timecourse aspect in the dataframe
            # FIXME: use xarray with timecourse dimension to handle the time dimension
            if k == 0:
                break

    # generate dataframes
    # TODO
    df_cells = pd.DataFrame(data_cells)
    df_points = pd.DataFrame(data_points)

    return df_cells, df_points






if __name__ == "__main__":
    console.rule(title="XDMF calculations", style="white")

    # FIXME: use the example of one of the zonation patterns with 10 timepoints
    xdmf_path = DATA_DIR / "spt" / "spt_zonation_patterns_new" / "simulation_pattern0.xdmf"
    xdmf_info: XDMFInformation = XDMFInformation.from_path(xdmf_path)
    console.print(xdmf_info)

    df_cells, df_points = create_mesh_dataframes(xdmf_path)
    console.rule(align="left", title="Cell data", style="white")
    console.print(df_cells)
    console.rule(align="left", title="point data", style="white")
    console.print(df_points)
