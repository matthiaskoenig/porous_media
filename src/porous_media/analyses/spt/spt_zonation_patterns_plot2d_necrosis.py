"""Visualization of liver simulations."""
from typing import Dict, List

import xarray as xr
from matplotlib import pyplot as plt

from porous_media import RESULTS_DIR
from porous_media.analyses.liver_variables import calculate_necrosis_fraction
from porous_media.console import console
from porous_media.data.xdmf_calculations import create_mesh_dataframes
from porous_media.visualization import plots


def plot_necrosis_over_time(
    xr_cells_list: List[xr.Dataset], labels: List[str], colors: List[str]
) -> None:
    """Plot necrosis over time."""
    console.rule(title="necrosis calculation", style="white")
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5, 5))

    # [1] necrosis fraction ~ time
    ax.set_xlabel("time [hr]", fontdict={"weight": "bold"})
    ax.set_ylabel("Necrosis [%]", fontdict={"weight": "bold"})

    for k, xr_cells in enumerate(xr_cells_list):
        necrosis_fraction = calculate_necrosis_fraction(xr_cells=xr_cells)

        ax.plot(
            # convert to hr and percent
            necrosis_fraction.time / 60 / 60,
            necrosis_fraction * 100,
            label=labels[k],
            linestyle="-",
            marker="o",
            color=colors[k],
            markeredgecolor="black",
            alpha=0.7,
        )
    ax.legend()
    plt.show()


if __name__ == "__main__":
    """Visualize 2D plots of the necrosis area."""

    # zonation analysis
    console.rule(title="XDMF calculations", style="white")
    # interpolated dataframe for zonation patterns
    xdmf_paths = [
        # RESULTS_DIR
        # / "spt_zonation_patterns_new"
        # / "10_28800.0"
        # / f"simulation_pattern{k}_interpolated.xdmf"
        RESULTS_DIR
        / "spt_zonation_patterns_new"
        / "100_28800.0"
        / f"simulation_pattern{k}_interpolated.xdmf"
        for k in range(5)
    ]

    # calculate all xarray Datasets
    xr_cells_list: List[xr.Dataset] = []
    for xdmf_path in xdmf_paths:
        xr_cells, xr_points = create_mesh_dataframes(xdmf_path)
        xr_cells_list.append(xr_cells)

    labels = [
        "Constant",
        "Linear increase",
        "Linear decrease",
        "Sharp pericentral",
        "Sharp periportal",
    ]
    colors = [
        "tab:blue",
        "tab:orange",
        "tab:green",
        "tab:red",
        "tab:purple",
    ]

    # calculate the necrosis area
    plot_necrosis_over_time(
        xr_cells_list=xr_cells_list,
        labels=labels,
        colors=colors,
    )
