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

    # DataFrame information
    df = simulation_conditions_df()

    console.rule(title="SPT timecourse", style="white")

    n_patterns = len(df.pattern_key.unique())
    n_cols = 4
    fig, axes = plt.subplots(nrows=n_patterns, ncols=n_cols, figsize=(n_cols*5, n_patterns*5), dpi=300)

    # [1] necrosis fraction ~ time
    for ax in axes.flatten():
        ax.set_xlabel("time [hr]", fontdict={"weight": "bold"})

    for k in range(n_patterns):
        axes[k, 0].set_ylabel("Substrate S [mM]", fontdict={"weight": "bold"})
        axes[k, 1].set_ylabel("Product P [mM]", fontdict={"weight": "bold"})
        axes[k, 2].set_ylabel("Toxic compound T [mM]", fontdict={"weight": "bold"})
        axes[k, 3].set_ylabel("Necrosis [%]", fontdict={"weight": "bold"})


    # TODO: identical y limits
    # TODO: boxplots with mean!
    # TODO: remove time axis
    # TODO: bigger axis labels
    # TODO: add title with pattern name
    # TODO: better colors, no orange

    ylim_maxs = {}
    for sim_id, xr_cells_raw in xr_cells_dict.items():

        info = df.loc[[sim_id]]
        color = info.color.values[0]
        pattern_key = info.pattern_key.values[0]

        # interpolate time
        xr_cells = xr_cells_raw.interp(time=times)

        kwargs = {
            "linestyle": "-",
            "marker": "o",
            "color": color,
            "markeredgecolor": "black",
        }

        for kax, sid in enumerate(["rr_(S)", "rr_(P)", "rr_(T)"]):


            x = xr_cells.time / 60 / 60,  # [s] -> [hr]
            y = xr_cells[sid].mean(dim="cell"),
            yerr = xr_cells[sid].std(dim="cell"),

            # update max
            if not sid in ylim_maxs:
                ylim_maxs[sid] = 0.0
            if (y + yerr).max()  > ylim_maxs[sid]:
                ylim_maxs[sid] = (y + yerr).max()

            ax = axes[pattern_key, kax]
            ax.errorbar(
                x=x,
                y=y,
                yerr=yerr,
                label=sim_id,
                **kwargs,
            )
            # ax.legend()

        necrosis_fraction = calculate_necrosis_fraction(xr_cells=xr_cells)
        axes[pattern_key, 3].plot(
            # convert to hr and percent
            necrosis_fraction.time / 60 / 60,  # [s] -> [hr]
            necrosis_fraction * 100,
            label=sim_id,
            **kwargs,
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

    from porous_media.analyses.spt.spt_information import pattern_names, pattern_colors, \
    simulation_conditions_df

    times: np.ndarray = np.linspace(start=0, stop=tend, num=51)
    plot_spt_over_time(
        xr_cells_dict=xr_cells_dict,
        times=times,
    )
