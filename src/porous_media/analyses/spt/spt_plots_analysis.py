"""Analysis of results."""
from pathlib import Path
from typing import Dict, List

import numpy as np
import xarray as xr
from matplotlib import pyplot as plt

from porous_media.analyses.liver_variables import calculate_necrosis_fraction
from porous_media.console import console
from porous_media.data.xdmf_calculations import mesh_datasets_from_xdmf

from porous_media.analyses.spt.spt_information import (
    pattern_idx2name,
    pattern_name2idx,
    pattern_order,
    boundary_flows,
    simulation_conditions_df,
)


def plot_positions(
    xr_cells_dict: Dict[str, xr.Dataset],
    xr_points_dict: Dict[str, xr.Dataset],
) -> None:
    """Plot positions."""
    console.rule(title="SPT positions", style="white")

    # DataFrame information
    df = simulation_conditions_df()
    n_patterns = len(pattern_order)
    n_cols = 6

    fig, axes = plt.subplots(nrows=n_patterns, ncols=n_cols, figsize=(n_cols*2.5, n_patterns*2.5),
                             dpi=300, layout="constrained")
    # [1] necrosis fraction ~ time

    for k_col in range(n_cols):
        axes[-1, k_col].set_xlabel("Position [-]", fontsize=9, fontdict={"weight": "bold"})

    for k_row in range(n_patterns):
        axes[k_row, 0].set_ylabel("Frequency [-]", fontsize=9, fontdict={"weight": "bold"})
        axes[k_row, 1].set_ylabel("Protein [-]", fontsize=9,
                                  fontdict={"weight": "bold"})
        axes[k_row, 2].set_ylabel("Mesh cell volume [nl]", fontsize=9, fontdict={"weight": "bold"})
        axes[k_row, 3].set_ylabel("Fluid volume [nl]", fontsize=9,
                                  fontdict={"weight": "bold"})
        axes[k_row, 4].set_ylabel("Solid volume [nl]", fontsize=9,
                                  fontdict={"weight": "bold"})
        axes[k_row, 5].set_ylabel("Pressure [?]", fontsize=9, fontdict={"weight": "bold"})

    ymax_vol = 0.0
    for k_row, pattern_name in enumerate(pattern_order):
        # only plot for the last pattern (full color)
        boundary_flow_key = len(boundary_flows) - 1
        pattern_key = pattern_name2idx[pattern_name]
        print(f"{pattern_name}, index={pattern_key}")
        df_sim = df[(df.pattern_key == pattern_key) & (df.boundary_flow_key == boundary_flow_key)]
        print(df_sim)

        sim_id = df_sim.index[0]
        color = df_sim.color.values[0]
        xr_cells_raw = xr_cells_dict[sim_id]
        xr_points_raw = xr_points_dict[sim_id]

        # Access last timepoint
        # FIXME: this is a bug should be the first time point, i.e. sel(time=0), isel(time=0)
        xr_cells: xr.Dataset = xr_cells_raw.isel(time=-1)
        xr_points: xr.Dataset = xr_points_raw.isel(time=-1)

        console.rule("cells")
        console.print(list(xr_cells.keys()))
        console.rule("points")
        console.print(list(xr_points.keys()))

        # Position histogram
        axes[k_row, 0].hist(
            xr_cells["rr_position"],
            color=color,
            edgecolor="black",
            alpha=0.7,
        )
        kwargs_scatter = {
            "marker": "o",
            "linestyle": "",
            "color": color,
            "markeredgecolor": "black",
            "alpha": 0.7,
        }
        # Protein scatter
        axes[k_row, 1].plot(
            xr_cells["rr_position"],
            xr_cells["rr_protein"],
            **kwargs_scatter,
        )
        # Volume scatter
        y = xr_cells["element_volume_TPM"]/1000/1E-9  # [m^3] -> nl
        ymax = y.max()
        if ymax > ymax_vol:
            ymax_vol = ymax
        axes[k_row, 2].plot(
            xr_cells["rr_position"],
            y,
            **kwargs_scatter,
        )
        # Fluid volume scatter
        y = xr_cells["rr_Vext"]/1E-9  # [l] -> [nl]
        ymax = y.max()
        if ymax > ymax_vol:
            ymax_vol = ymax
        axes[k_row, 3].plot(
            xr_cells["rr_position"],
            y,
            **kwargs_scatter,
        )
        # Solid volume scatter
        y = xr_cells["rr_Vli"]/1E-9  # [l] -> [nl]
        ymax = y.max()
        if ymax > ymax_vol:
            ymax_vol = ymax
        axes[k_row, 4].plot(
            xr_cells["rr_position"],
            xr_cells["rr_Vli"]/1E-9,
            **kwargs_scatter,
        )
        # Pressure scatter
        axes[k_row, 5].plot(
            xr_cells["rr_position"],
            xr_cells["pressure"],
            **kwargs_scatter,
        )
        console.print(xr_cells["pressure"].values)

    for k_row, pattern_name in enumerate(pattern_order):
        axes[k_row, 0].set_title(pattern_name, fontsize=15, fontweight="bold")
        for k_col in range(n_cols):
            axes[k_row, k_col].set_xlim([0, 1])
        for k_col in range(n_cols-1):
            axes[k_row, k_col].set_ylim(bottom=0)

        for k_col in [2, 3, 4]:
            axes[k_row, k_col].set_ylim(top=1.05*ymax_vol)

    plt.show()


def plot_spt_over_time(
    xr_cells_dict: Dict[str, xr.Dataset],
    times: np.ndarray,
) -> None:
    """Plot SPT over time."""
    console.rule(title="SPT timecourse", style="white")

    # DataFrame information
    df = simulation_conditions_df()
    n_patterns = len(pattern_order)
    n_cols = 6

    fig, axes = plt.subplots(nrows=n_patterns, ncols=n_cols, figsize=(n_cols*2.5, n_patterns*2.5),
                             dpi=300, layout="constrained")
    # [1] necrosis fraction ~ time
    for ax in axes[-1, :].flatten():
        ax.set_xlabel("time [hr]", fontsize=9, fontdict={"weight": "bold"})

    for k_row in range(n_patterns):
        axes[k_row, 0].set_ylabel("Substrate ext [mM]", fontsize=9, fontdict={"weight": "bold"})
        axes[k_row, 1].set_ylabel("Product ext [mM]", fontsize=9, fontdict={"weight": "bold"})
        axes[k_row, 2].set_ylabel("Substrate [mM]", fontsize=9, fontdict={"weight": "bold"})
        axes[k_row, 3].set_ylabel("Product [mM]", fontsize=9, fontdict={"weight": "bold"})
        axes[k_row, 4].set_ylabel("Toxic [mM]", fontsize=9, fontdict={"weight": "bold"})
        axes[k_row, 5].set_ylabel("Necrosis [%]", fontsize=9, fontdict={"weight": "bold"})
        axes[k_row, 5].set_ylim([0, 100*1.05])

    # TODO: boxplots with mean!

    ylim_maxs = {}
    for k_row, pattern_name in enumerate(pattern_order):
        for boundary_flow_key in range(len(boundary_flows)):
            pattern_key = pattern_name2idx[pattern_name]
            df_sim = df[(df.pattern_key == pattern_key) & (df.boundary_flow_key == boundary_flow_key)]
            sim_id = df_sim.index[0]
            color = df_sim.color.values[0]
            xr_cells_raw = xr_cells_dict[sim_id]

            # interpolate time
            xr_cells = xr_cells_raw.interp(time=times)

            kwargs = {
                "linestyle": "-",
                "marker": "o",
                "color": color,
                "markeredgecolor": "black",
            }

            for k_col, sid in enumerate(["rr_(S_ext)", "rr_(P_ext)", "rr_(S)", "rr_(P)", "rr_(T)"]):

                x = xr_cells.time / 60 / 60  # [s] -> [hr]
                y = xr_cells[sid].mean(dim="cell")
                yerr = xr_cells[sid].std(dim="cell")

                # update max
                if not sid in ylim_maxs:
                    ylim_maxs[sid] = 0.0
                if (y + yerr).max() > ylim_maxs[sid]:
                    ylim_maxs[sid] = (y + yerr).max()

                ax = axes[k_row, k_col]
                ax.errorbar(
                    x=x,
                    y=y,
                    yerr=yerr,
                    label=sim_id,
                    **kwargs,
                )
                # ax.legend()

            necrosis_fraction = calculate_necrosis_fraction(xr_cells=xr_cells)
            axes[k_row, 5].plot(
                # convert to hr and percent
                necrosis_fraction.time / 60 / 60,  # [s] -> [hr]
                necrosis_fraction * 100,
                label=sim_id,
                **kwargs,
            )

    for k_row, pattern_name in enumerate(pattern_order):
        axes[k_row, 0].set_title(pattern_name, fontsize=15, fontweight="bold")
        for kax, sid in enumerate(["rr_(S_ext)", "rr_(P_ext)", "rr_(S)", "rr_(P)", "rr_(T)"]):
            axes[k_row, kax].set_ylim([0, 1.05*ylim_maxs[sid]])

    plt.show()


def plot_spt_over_position(
    xr_cells_dict: Dict[str, xr.Dataset],
) -> None:
    """Plot SPT over position."""
    console.rule(title="SPT position", style="white")
    df = simulation_conditions_df()
    n_patterns = len(pattern_order)
    n_cols = 6
    fig, axes = plt.subplots(nrows=n_patterns, ncols=n_cols, figsize=(n_cols*2.5, n_patterns*2.5),
                             dpi=300, layout="constrained")

    for ax in axes[-1, :].flatten():
        ax.set_xlabel("Position [-]", fontsize=9, fontdict={"weight": "bold"})

    for k_row in range(n_patterns):
        axes[k_row, 0].set_ylabel("Substrate ext [mM]", fontsize=9, fontdict={"weight": "bold"})
        axes[k_row, 1].set_ylabel("Product ext [mM]", fontsize=9, fontdict={"weight": "bold"})
        axes[k_row, 2].set_ylabel("Substrate [mM]", fontsize=9, fontdict={"weight": "bold"})
        axes[k_row, 3].set_ylabel("Product [mM]", fontsize=9, fontdict={"weight": "bold"})
        axes[k_row, 4].set_ylabel("Toxic [mM]", fontsize=9, fontdict={"weight": "bold"})
        axes[k_row, 5].set_ylabel("Necrosis [-]", fontsize=9, fontdict={"weight": "bold"})

    # TODO: boxplots with mean!

    ylim_maxs = {}
    for k_row, pattern_name in enumerate(pattern_order):
        for boundary_flow_key in range(len(boundary_flows)):
            pattern_key = pattern_name2idx[pattern_name]
            df_sim = df[(df.pattern_key == pattern_key) & (df.boundary_flow_key == boundary_flow_key)]
            sim_id = df_sim.index[0]
            color = df_sim.color[0]
            xr_cells_raw = xr_cells_dict[sim_id]

            # interpolate time (only last timepoint)
            xr_cells = xr_cells_raw.interp(time=10*60*60)  # 10 hr

            kwargs = {
                "linestyle": "",
                "marker": "o",
                "color": color,
                "markeredgecolor": "black",
            }

            for k_col, sid in enumerate(["rr_(S_ext)", "rr_(P_ext)", "rr_(S)", "rr_(P)", "rr_(T)", "rr_necrosis"]):
                x = xr_cells.rr_position
                y = xr_cells[sid]

                # update max
                if not sid in ylim_maxs:
                    ylim_maxs[sid] = 0.0
                if y.max() > ylim_maxs[sid]:
                    ylim_maxs[sid] = y.max()

                ax = axes[k_row, k_col]
                ax.plot(
                    x,
                    y,
                    label=sim_id,
                    **kwargs,
                )
                # ax.legend()

    for k_row, pattern_name in enumerate(pattern_order):
        axes[k_row, 0].set_title(pattern_name, fontsize=15, fontweight="bold")
        for k_col, sid in enumerate(["rr_(S_ext)", "rr_(P_ext)", "rr_(S)", "rr_(P)", "rr_(T)", "rr_necrosis"]):
            axes[k_row, k_col].set_ylim([0, 1.05*ylim_maxs[sid]])

    plt.show()


if __name__ == "__main__":
    """Analysis plots of the TMP simulations."""

    # TODO: add protein to visualization
    # TODO: fix the zonation axes to 0, 1
    # TODO: create boxplots


    # zonation analysis
    console.rule(title="Dataset calculation", style="white")

    # XDMF
    xdmf_dir = Path("/home/mkoenig/git/porous_media/data/spt/2023-12-06/xdmf")
    xdmf_paths = [f for f in xdmf_dir.glob("*.xdmf")]

    # Load xarray datasets
    xr_cells_dict: Dict[str, xr.Dataset] = {}
    xr_points_dict: Dict[str, xr.Dataset] = {}
    tend: float = np.Inf
    for xdmf_path in xdmf_paths:
        xr_cells, xr_points = mesh_datasets_from_xdmf(xdmf_path)

        # get sim id and store data
        sim_id = xdmf_path.stem
        xr_cells_dict[sim_id] = xr_cells
        xr_points_dict[sim_id] = xr_points

        tend_sim = xr_cells.time[-1]
        print(tend_sim)
        if tend_sim < tend:
            tend = tend_sim

    # figure out end time
    times: np.ndarray = np.linspace(start=0, stop=tend, num=51)

    plot_positions(
        xr_cells_dict=xr_cells_dict,
        xr_points_dict=xr_points_dict,
    )

    # plot_spt_over_time(
    #     xr_cells_dict=xr_cells_dict,
    #     times=times,
    # )
    #
    # plot_spt_over_position(
    #     xr_cells_dict=xr_cells_dict,
    # )
