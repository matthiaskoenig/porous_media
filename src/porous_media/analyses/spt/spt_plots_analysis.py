"""Analysis of results."""
from pathlib import Path
from typing import Dict, List

import numpy as np
import xarray as xr
from matplotlib import pyplot as plt

from porous_media.analyses.liver_variables import calculate_necrosis_fraction
from porous_media.console import console
from porous_media.data.xdmf_calculations import mesh_datasets_from_xdmf

def plot_spt_over_time(
    xr_cells_dict: Dict[str, xr.Dataset],
    times: np.ndarray,
) -> None:
    """Plot SPT over time."""
    console.rule(title="necrosis calculation", style="white")
    fig, axes = plt.subplots(nrows=1, ncols=4, figsize=(20, 5))

    # [1] necrosis fraction ~ time
    for ax in axes:
        ax.set_xlabel("time [hr]", fontdict={"weight": "bold"})

    axes[0].set_ylabel("Substrate S [mM]", fontdict={"weight": "bold"})
    axes[1].set_ylabel("Product P [mM]", fontdict={"weight": "bold"})
    axes[2].set_ylabel("Toxic compound T [mM]", fontdict={"weight": "bold"})
    axes[3].set_ylabel("Necrosis [%]", fontdict={"weight": "bold"})

    for sim_id, xr_cells_raw in xr_cells_dict.items():

        # interpolate time
        xr_cells = xr_cells_raw.interp(time=times)

        for kax, sid in enumerate(["rr_(S)", "rr_(P)", "rr_(T)"]):
            ax = axes[kax]
            ax.errorbar(
                # convert to hr and percent
                x=xr_cells.time / 60 / 60,  # [s] -> [hr]
                y=xr_cells[sid].mean(dim="cell"),
                yerr=xr_cells[sid].std(dim="cell"),
                label=sim_id,
                linestyle="-",
                marker="o",
                # color=colors[k],
                # markeredgecolor="black",
                alpha=0.7,
            )
            ax.legend()

        necrosis_fraction = calculate_necrosis_fraction(xr_cells=xr_cells)
        axes[3].plot(
            # convert to hr and percent
            necrosis_fraction.time / 60 / 60, # [s] -> [hr]
            necrosis_fraction * 100,
            label=sim_id,
            linestyle="-",
            marker="o",
            # color=colors[k],
            # markeredgecolor="black",
            alpha=0.7,
        )

    plt.show()


if __name__ == "__main__":
    """Analysis plots of the TMP simulations."""

    # zonation analysis
    console.rule(title="Dataset calculation", style="white")

    # interpolated dataframe for zonation patterns
    xdmf_dir = Path("/home/mkoenig/git/porous_media/data/spt/2023-12-05")
    xdmf_paths = [f for f in xdmf_dir.glob("*.xdmf")]

    # xdmf_paths = [
    #     xdmf_dir / "sim001.xdmf",
    #     xdmf_dir / "sim0010.xdmf",
    #     xdmf_dir / "sim0011.xdmf",
    #     xdmf_dir / "sim0012.xdmf",
    # ]
    # load all xarray Datasets
    xr_cells_dict: Dict[str, xr.Dataset] = {}
    tend: float = np.Inf
    for xdmf_path in xdmf_paths:
        xr_cells, xr_points = mesh_datasets_from_xdmf(xdmf_path)
        sim_id = xdmf_path.stem
        xr_cells_dict[sim_id] = xr_cells
        tend_sim = xr_cells.time[-1]
        print(tend_sim)
        if tend_sim < tend:
            tend = tend_sim


    # figure out end time

    from porous_media.analyses.spt.spt_information import pattern_names, pattern_colors

    times: np.ndarray = np.linspace(start=0, stop=tend, num=51)
    plot_spt_over_time(
        xr_cells_dict=xr_cells_dict,
        times=times,
    )
