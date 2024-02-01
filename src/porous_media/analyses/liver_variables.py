"""Helper module for calculating variables from the simulation results.

Here reused variables should be defined.
"""

import xarray as xr

from porous_media.console import console


def calculate_necrosis_fraction(xr_cells: xr.Dataset) -> xr.Dataset:
    """Calculate the necrosis fraction from a given geometry.

    The necrosis and cell volumes are used to calculate the fraction.

    returns necrosis_fraction in [0, 1]
    """

    if not hasattr(xr_cells, "rr_necrosis"):
        raise ValueError(f"No attribute/variable 'necrosis' in cell data: {xr_cells}")

    necrosis = xr_cells.rr_necrosis
    # replace necrosis values > 0 with 1.0 (partial necrosis is necrosis)
    # # FIXME
    # console.print("before:")
    # console.print(necrosis)
    # necrosis = necrosis.where(necrosis > 1E-3, 1.0)
    # console.print("after:")

    if "element_volume_point_TPM" in xr_cells:
        cell_volumes = xr_cells.element_volume_point_TPM
    else:
        # FIXME: bugfix
        cell_volumes = xr_cells.element_volume_TPM
    if necrosis.min() < 0.0:
        raise ValueError(f"'necrosis' must be => 0.0, but minimum is {necrosis.min()}")
    if necrosis.max() > 1.0:
        raise ValueError(f"'necrosis' must be <= 1.0, but maximum is {necrosis.max()}")

    # calculate necrosis fraction over time
    necrosis_fraction = (necrosis * cell_volumes).sum(dim="cell") / cell_volumes.sum(
        dim="cell"
    )

    return necrosis_fraction
