"""Helper module for calculating variables from the simulation results.

Here reused variables should be defined.
"""
from pathlib import Path
from typing import List

import xarray as xr

from porous_media.console import console
from porous_media.data.xdmf_calculations import create_mesh_dataframes
from porous_media.data.xdmf_tools import XDMFInfo


def calculate_cell_geometries(xr_cells: xr.Dataset) -> None:
    """Calculate the geometrical parameters of the cells in the mesh.

    This adds/could add additional variables to the mesh for analysis.
    """

    # FIXME: not implemented



    return [1.0]


def calculate_total_volume() -> float:
    """Calculate total volume of geometry."""

    # FIXME: implement; this is also time dependent
    # see
    # https://stackoverflow.com/questions/61638966/find-volume-of-object-given-a-triangular-mesh
    # V = np.sum(mesh.cell_volumes)
    # print("Volume =", V)
    # Volume = 5.6118446
    raise NotImplementedError
    return 1.0


def calculate_necrosis_fraction(xr_cells: xr.Dataset) -> xr.Dataset:
    """Calculate the necrosis fraction from a given geometry.

    returns necrosis_fraction in [0, 1], to scale to percentage volume
    """

    if not hasattr(xr_cells, "rr_necrosis"):
        raise ValueError(f"No attribute/variable 'necrosis' in cell data: {xr_cells}")

    necrosis = xr_cells.rr_necrosis
    if necrosis.min() < 0.0:
        raise ValueError(f"'necrosis' must be => 0.0, but minimum is {necrosis.min()}")
    if necrosis.max() > 1.0:
        raise ValueError(f"'necrosis' must be <= 1.0, but maximum is {necrosis.max()}")

    # calculate necrosis fraction (sum/count)
    necrosis_fraction = necrosis.sum(dim="cell") / necrosis.count(dim="cell")
    # FIXME: calculate and add the cell volumes for proper normalization, i.e., the
    # different cells have different volumes.

    return necrosis_fraction


if __name__ == "__main__":
    # simple example for adding variables (load interpolated xdmf)
    from porous_media import RESULTS_DIR

    xdmf_path: Path = RESULTS_DIR / "spt_zonation_patterns_219/10_28800.0/simulation_pattern1_interpolated.xdmf"
    xdmf_info = XDMFInfo.from_path(xdmf_path)
    console.print(xdmf_info)

    console.rule(style="white")
    xr_cells: xr.Dataset
    xr_points: xr.Dataset
    xr_cells, xr_points = create_mesh_dataframes(xdmf_path)
    console.print(xr_cells)


