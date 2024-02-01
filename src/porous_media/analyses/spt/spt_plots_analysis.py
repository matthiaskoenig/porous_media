"""Analysis of results."""

from pathlib import Path
from typing import Dict, List

import numpy as np
import xarray as xr
from matplotlib import pyplot as plt

from porous_media.analyses.liver_variables import calculate_necrosis_fraction
from porous_media.analyses.spt.spt_information import (
    boundary_flows,
    pattern_name2idx,
    pattern_order,
    simulation_conditions_df,
)
from porous_media.console import console
from porous_media.data.xdmf_calculations import mesh_datasets_from_xdmf


label_kwargs = {
    "fontsize": 12,
    "fontweight": "bold",
}


def plot_positions(
    results_dir: Path,
    xr_cells_dict: Dict[str, xr.Dataset],
    xr_points_dict: Dict[str, xr.Dataset],
) -> None:
    """Plot positions."""
    console.rule(title="SPT positions", style="white")

    # DataFrame information
    df = simulation_conditions_df()
    n_patterns = len(pattern_order)
    n_cols = 6

    fig, axes = plt.subplots(
        nrows=n_patterns,
        ncols=n_cols,
        figsize=(n_cols * 2.7, n_patterns * 2.5),
        dpi=300,
        layout="constrained",
    )
    # [1] necrosis fraction ~ time

    for k_col in range(n_cols):
        axes[-1, k_col].set_xlabel("Position [-]", **label_kwargs)

    for k_row in range(n_patterns):
        axes[k_row, 0].set_ylabel("Frequency [-]", **label_kwargs)
        axes[k_row, 1].set_ylabel("Protein [-]", **label_kwargs)
        axes[k_row, 2].set_ylabel("Mesh cell volume [nl]", **label_kwargs)
        axes[k_row, 3].set_ylabel("Fluid volume [nl]", **label_kwargs)
        axes[k_row, 4].set_ylabel("Solid volume [nl]", **label_kwargs)
        axes[k_row, 5].set_ylabel("Pressure [?]", **label_kwargs)

    ymax_vol = 0.0
    for k_row, pattern_name in enumerate(pattern_order):
        # only plot for the last pattern (full color)
        boundary_flow_key = len(boundary_flows) - 1
        pattern_key = pattern_name2idx[pattern_name]
        print(f"{pattern_name}, index={pattern_key}")
        df_sim = df[
            (df.pattern_key == pattern_key)
            & (df.boundary_flow_key == boundary_flow_key)
        ]
        print(df_sim)

        sim_id = df_sim.index[0]
        color = df_sim.color.values[0]
        xr_cells_raw = xr_cells_dict[sim_id]
        xr_points_raw = xr_points_dict[sim_id]

        # Access last timepoint
        # FIXME: should be first point with all the data
        xr_cells: xr.Dataset = xr_cells_raw.isel(time=-1)
        _: xr.Dataset = xr_points_raw.isel(time=-1)

        console.rule("cells")
        console.print(list(xr_cells.keys()))

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
            # "markeredgewidth": 0.5,
            "alpha": 0.7,
        }
        # Protein scatter
        axes[k_row, 1].plot(
            xr_cells["rr_position"],
            xr_cells["rr_protein"],
            **kwargs_scatter,
        )
        # Volume scatter
        y = xr_cells["element_volume_point_TPM"] / 1e-12  # [m^3] -> [nl]
        ymax = y.max()
        if ymax > ymax_vol:
            ymax_vol = ymax
        axes[k_row, 2].plot(
            xr_cells["rr_position"],
            y,
            **kwargs_scatter,
        )

        # Fluid volume scatter
        y = xr_cells["rr_Vext"] / 1e-9  # [l] -> [nl]
        ymax = y.max()
        if ymax > ymax_vol:
            ymax_vol = ymax
        axes[k_row, 3].plot(
            xr_cells["rr_position"],
            y,
            **kwargs_scatter,
        )
        # Solid volume scatter
        y = xr_cells["rr_Vli"] / 1e-9  # [l] -> [nl]
        ymax = y.max()
        if ymax > ymax_vol:
            ymax_vol = ymax
        axes[k_row, 4].plot(
            xr_cells["rr_position"],
            xr_cells["rr_Vli"] / 1e-9,
            **kwargs_scatter,
        )
        # Pressure scatter
        axes[k_row, 5].plot(
            xr_cells["rr_position"],
            xr_cells["pressure"],
            **kwargs_scatter,
        )

    for k_row, _ in enumerate(pattern_order):
        # axes[k_row, 0].set_title(pattern_name, fontsize=15, fontweight="bold")

        for k_col in range(n_cols):
            axes[k_row, k_col].set_xlim([0, 1])
        for k_col in range(n_cols - 1):
            axes[k_row, k_col].set_ylim(bottom=0)

        for k_col in [3, 4]:
            axes[k_row, k_col].set_ylim(bottom=-0.05 * ymax_vol, top=1.05 * ymax_vol)

    for ax in axes[:-1, :].flatten():
        ax.set_xticks([0, 0.5, 1.0])
        ax.set_xticklabels([])
    for ax in axes[-1, :].flatten():
        ax.xaxis.set_ticks([0, 0.5, 1], labels=["PP", "0.5", "PV"])

    plt.show()
    fig.savefig(results_dir / "position_tests.png", bbox_inches="tight")


def plot_spt_over_time(
    results_dir: Path,
    xr_cells_dict: Dict[str, xr.Dataset],
    times: np.ndarray,
) -> None:
    """Plot SPT over time."""
    console.rule(title="SPT timecourse", style="white")

    # DataFrame information
    df = simulation_conditions_df()
    n_patterns = len(pattern_order)
    n_cols = 5

    fig, axes = plt.subplots(
        nrows=n_patterns,
        ncols=n_cols,
        figsize=(n_cols * 2.7, n_patterns * 2.5),
        dpi=300,
        layout="constrained",
    )
    # [1] necrosis fraction ~ time
    for ax in axes[-1, :].flatten():
        ax.set_xlabel("Time [hr]", **label_kwargs)

    for k_row in range(n_patterns):
        axes[k_row, 0].set_ylabel("Substrate [mM]", **label_kwargs)
        axes[k_row, 1].set_ylabel("Product [mM]", **label_kwargs)
        axes[k_row, 2].set_ylabel("Protein [-]", **label_kwargs)
        axes[k_row, 3].set_ylabel("Toxic compound [mM]", **label_kwargs)
        axes[k_row, 4].set_ylabel("Necrosis [%]", **label_kwargs)
        axes[k_row, 4].set_ylim([0, 100 * 1.05])

    ylim_maxs = {}
    for k_row, pattern_name in enumerate(pattern_order):
        for boundary_flow_key in range(len(boundary_flows)):
            pattern_key = pattern_name2idx[pattern_name]
            df_sim = df[
                (df.pattern_key == pattern_key)
                & (df.boundary_flow_key == boundary_flow_key)
            ]
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
                # "markeredgewidth": 0.5,
            }

            for k_col, sid in enumerate(
                ["rr_(S_ext)", "rr_(P_ext)", "rr_protein", "rr_(T)"]
            ):
                x = xr_cells.time / 60 / 60  # [s] -> [hr]
                y = xr_cells[sid].mean(dim="cell")
                yerr = xr_cells[sid].std(dim="cell")

                # update max
                if sid not in ylim_maxs:
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
            axes[k_row, 4].plot(
                # convert to hr and percent
                necrosis_fraction.time / 60 / 60,  # [s] -> [hr]
                necrosis_fraction * 100,
                label=sim_id,
                **kwargs,
            )

    for k_row, _ in enumerate(pattern_order):
        # axes[k_row, 2].set_title(pattern_name, fontsize=20, fontweight="bold")

        for kax, sid in enumerate(["rr_(S_ext)", "rr_(P_ext)", "rr_protein", "rr_(T)"]):
            axes[k_row, kax].set_ylim([-0.05 * ylim_maxs[sid], 1.05 * ylim_maxs[sid]])

    for ax in axes[:-1, :].flatten():
        ax.set_xticks([0, 10, 20])
        ax.set_xticklabels([])
    for ax in axes[-1, :].flatten():
        ax.xaxis.set_ticks([0, 10, 20], labels=["0", "10", "20"])

    plt.show()
    fig.savefig(results_dir / "time_dependency.png", bbox_inches="tight")


def plot_spt_over_position(
    results_dir: Path,
    xr_cells_dict: Dict[str, xr.Dataset],
) -> None:
    """Plot SPT over position."""
    console.rule(title="SPT position", style="white")
    df = simulation_conditions_df()
    n_patterns = len(pattern_order)
    n_cols = 5
    fig, axes = plt.subplots(
        nrows=n_patterns,
        ncols=n_cols,
        figsize=(n_cols * 2.7, n_patterns * 2.5),
        dpi=300,
        layout="constrained",
    )

    for ax in axes[-1, :].flatten():
        ax.set_xlabel("Position [-]", **label_kwargs)

    for k_row in range(n_patterns):
        axes[k_row, 0].set_ylabel("Substrate [mM]", **label_kwargs)
        axes[k_row, 1].set_ylabel("Product [mM]", **label_kwargs)
        axes[k_row, 2].set_ylabel("Protein [-]", **label_kwargs)
        axes[k_row, 3].set_ylabel("Toxic compound [mM]", **label_kwargs)
        axes[k_row, 4].set_ylabel("Necrosis [-]", **label_kwargs)

    ylim_maxs = {}
    for k_row, pattern_name in enumerate(pattern_order):
        for boundary_flow_key in range(len(boundary_flows)):
            pattern_key = pattern_name2idx[pattern_name]
            df_sim = df[
                (df.pattern_key == pattern_key)
                & (df.boundary_flow_key == boundary_flow_key)
            ]
            sim_id = df_sim.index[0]
            color = df_sim.color[0]

            xr_cells_raw = xr_cells_dict[sim_id]

            # interpolate time (only last timepoint)
            xr_cells = xr_cells_raw.interp(time=10 * 60 * 60)  # 10 hr

            kwargs = {
                "linestyle": "",
                "marker": "o",
                "color": color,
                "markeredgecolor": "black",
            }

            for k_col, sid in enumerate(
                ["rr_(S_ext)", "rr_(P_ext)", "rr_protein", "rr_(T)", "rr_necrosis"]
            ):
                x = xr_cells.rr_position
                y = xr_cells[sid]

                # filter necrosis values != 0.0 or 1.0 (partial necrosis due to point averaging)
                if sid == "rr_necrosis":
                    x = x.where((y == 0.0) | (y == 1.0))
                    y = y.where((y == 0.0) | (y == 1.0))

                # update max
                if sid not in ylim_maxs:
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

    for k_row, _ in enumerate(pattern_order):
        # axes[k_row, 2].set_title(pattern_name, fontsize=20, fontweight="bold")
        for k_col, sid in enumerate(
            ["rr_(S_ext)", "rr_(P_ext)", "rr_protein", "rr_(T)", "rr_necrosis"]
        ):
            axes[k_row, k_col].set_ylim([-0.05 * ylim_maxs[sid], 1.05 * ylim_maxs[sid]])

    for ax in axes[:-1, :].flatten():
        ax.set_xticks([0, 0.5, 1.0])
        ax.set_xticklabels([])
    for ax in axes[-1, :].flatten():
        ax.xaxis.set_ticks([0, 0.5, 1], labels=["PP", "", "PV"])

    plt.show()
    fig.savefig(results_dir / "position_dependency.png", bbox_inches="tight")


if __name__ == "__main__":
    """Analysis plots of the SPT simulations."""

    # date = "2023-12-13"
    date = "2023-12-19"
    console.rule(title=f"SPT analysis: {date}", style="white")

    # XDMF
    xdmf_dir = Path(f"/home/mkoenig/git/porous_media/data/spt/{date}/xdmf")
    xdmf_paths = sorted([f for f in xdmf_dir.glob("*.xdmf")])

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
        if tend_sim < tend:
            tend = tend_sim

    # figure out end time
    times: np.ndarray = np.linspace(start=0, stop=tend, num=51)

    from porous_media import RESULTS_DIR

    results_dir = RESULTS_DIR / "spt" / date
    results_dir.mkdir(exist_ok=True, parents=True)

    plot_positions(
        results_dir=results_dir,
        xr_cells_dict=xr_cells_dict,
        xr_points_dict=xr_points_dict,
    )

    plot_spt_over_time(
        results_dir=results_dir,
        xr_cells_dict=xr_cells_dict,
        times=times,
    )

    plot_spt_over_position(
        results_dir=results_dir,
        xr_cells_dict=xr_cells_dict,
    )
