"""
Visualization of liver
"""
from typing import Dict

import xarray as xr

from matplotlib import pyplot as plt

from porous_media import RESULTS_DIR
from porous_media.console import console
from porous_media.data.xdmf_calculations import create_mesh_dataframes


def plot_necrosis_over_time(xr_cells_dict: Dict[str, xr.Dataset]):
    """Plot necrosis over time."""
    console.rule(title="necrosis calculation", style="white")
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7, 7))

    # [1] necrosis fraction ~ time
    ax.set_xlabel("time [hr]")
    ax.set_ylabel("necrosis fraction [%]")

    for label, xr_cells in xr_cells_dict.items():
        # calculate necrosis fraction (sum/count)
        # FIXME: calculate and add the cell volumes for proper normalization
        necrosis = xr_cells.rr_necrosis
        necrosis_fraction = necrosis.sum(dim="cell")/necrosis.count(dim="cell")
        console.print(f"{necrosis_fraction=}")

        ax.plot(
            # convert to hr and percent
            necrosis_fraction.time/60/60, necrosis_fraction * 100,
            label=label,
            linestyle="-",
            marker="o",
            markeredgecolor="black",
            alpha=0.7
        )
    ax.legend()
    plt.show()


if __name__ == "__main__":
    # zonation analysis
    console.rule(title="XDMF calculations", style="white")
    # interpolated dataframe for zonation patterns
    xdmf_paths = [
        RESULTS_DIR / "spt_zonation_patterns_new" / "10_28800.0" / f"simulation_pattern{k}_interpolated.xdmf"
        # RESULTS_DIR / "spt_zonation_patterns_new" / "100_28800.0" / f"simulation_pattern{k}_interpolated.xdmf"
        for k in range(5)]
    labels = [
        "Constant",
        "Linear increase",
        "Linear decrease",
        "Sharp pericentral",
        "Sharp periportal",
    ]

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