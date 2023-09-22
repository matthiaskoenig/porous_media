"""
Visualization of liver
"""
from typing import Dict, List

import xarray as xr
from porous_media.visualization import plots
from matplotlib import pyplot as plt

from porous_media import RESULTS_DIR
from porous_media.console import console
from porous_media.data.xdmf_calculations import create_mesh_dataframes


def plot_necrosis_over_time(
    xr_cells_list: List[xr.Dataset], labels: List[str], colors: List[str]
):
    """Plot necrosis over time."""
    console.rule(title="necrosis calculation", style="white")
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5, 5))
    fig.subplots_adjust(wspace=0.3, hspace=0.3)

    # [1] necrosis fraction ~ time
    ax.set_xlabel("time [hr]", fontdict={"weight": "bold"})
    ax.set_ylabel("Necrosis [%]", fontdict={"weight": "bold"})

    for k, xr_cells in enumerate(xr_cells_list):
        # calculate necrosis fraction (sum/count)
        # FIXME: calculate and add the cell volumes for proper normalization
        necrosis = xr_cells.rr_necrosis
        necrosis_fraction = necrosis.sum(dim="cell") / necrosis.count(dim="cell")
        # console.print(f"{necrosis_fraction=}")

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
    fig.savefig("flow_dependency.png", dpi=300)
    plt.show()



def plot_T_vs_position(
    xr_cells_list: List[xr.Dataset], labels: List[str], colors: List[str]
):
    """Plot necrosis vs position."""
    console.rule(title="necrosis calculation", style="white")
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7, 7))
    fig.adjust_subplots(wspace=0.3, hspace=0.3)

    # [1] necrosis fraction ~ time
    ax.set_xlabel("position [-]")
    ax.set_ylabel("toxic compound [mM]")

    for k, xr_cells in enumerate(xr_cells_list):
        ax.plot(
            # convert to hr and percent
            xr_cells.rr_position,
            xr_cells["rr_(T)"],
            # label=labels[k],
            linestyle="None",
            marker="o",
            color=colors[k],
            markeredgecolor="black",
            alpha=0.7,
        )
    plt.show()


if __name__ == "__main__":
    # zonation analysis
    console.rule(title="XDMF calculations", style="white")
    # interpolated dataframe for zonation patterns
    xdmf_paths = [
        # RESULTS_DIR
        # / "spt_zonation_patterns_new"
        # / "10_28800.0"
        # / f"simulation_pattern{k}_interpolated.xdmf"
        RESULTS_DIR / "spt_substrate_scan_219" / "100_28800.0" / f"spt_{k}_interpolated.xdmf"
        for k in range(21, 26)
    ]

    labels = [
        "flow S(1)",
        "flow S(2)",
        "flow S(3)",
        "flow S(4)",
        "flow S(5)",
    ]
    colors = [
        "tab:blue",
        "tab:orange",
        "tab:green",
        "tab:red",
        "tab:purple",
    ]

    # calculate all xarray Datasets
    xr_cells_list: List[xr.Dataset] = []
    for k, xdmf_path in enumerate(xdmf_paths):
        xr_cells, xr_points = create_mesh_dataframes(xdmf_path)
        xr_cells_list.append(xr_cells)

    # calculate the necrosis area
    plot_necrosis_over_time(
        xr_cells_list=xr_cells_list,
        labels=labels,
        colors=colors,
    )

    plot_T_vs_position(
        xr_cells_list=xr_cells_list,
        labels=labels,
        colors=colors,
    )
