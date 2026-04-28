"""Analysis of convergence."""

from pathlib import Path

import matplotlib
import numpy as np
import xarray as xr
from matplotlib import pyplot as plt

from porous_media.console import console
from porous_media.data.xdmf_calculations import mesh_datasets_from_xdmf


label_kwargs = {
    "fontsize": 12,
    "fontweight": "bold",
}

convergence_resolutions = [
    # "00005",
    "000025",
    "000015",
    "0000125",
    "00001",
    "00000625",
]
convergence_labels = {
    # "00005": "5.000",
    "000025": "2.500",
    "000015": "1.500",
    "0000125": "1.250",
    "00001": "1.000",
    "00000625": "0.625",
}


def f_convergence_colors() -> dict[str, str]:
    """Calculate colors."""
    colors = {}
    cmap = matplotlib.colormaps.get_cmap("viridis")
    n_resolutions = len(convergence_resolutions)
    for kres, resolution in enumerate(convergence_resolutions):
        color_rgba = cmap(kres / (n_resolutions - 1))
        colors[resolution] = matplotlib.colors.to_hex(color_rgba, keep_alpha=True)
    return colors


convergence_colors: dict[str, str] = f_convergence_colors()

kwargs = {
    "linestyle": "-",
    "marker": "None",  # "o",
    "markeredgecolor": "black",
    "markeredgewidth": 0.5,
    "markersize": 3,
}


def plot_convergence_over_time(
    results_dir: Path,
    xr_cells_dict: dict[str, xr.Dataset],
    times: np.ndarray,
) -> None:
    """Plot convergence over time."""
    console.rule(title="timecourse", style="white")

    # DataFrame information
    n_cols = len(convergence_resolutions)
    fig, axes = plt.subplots(
        ncols=n_cols,
        figsize=(n_cols * 2.7, 2.5),
        dpi=300,
        layout="constrained",
    )
    # [1] necrosis fraction ~ time
    for ax in axes:
        ax.set_xlabel("Time [hr]", **label_kwargs)

    axes[0].set_ylabel("Substrate [mM]", **label_kwargs)
    axes[1].set_ylabel("Product [mM]", **label_kwargs)
    axes[2].set_ylabel("Protein [-]", **label_kwargs)
    axes[3].set_ylabel("Toxic compound [mM]", **label_kwargs)
    axes[4].set_ylabel("Necrosis [%]", **label_kwargs)
    # axes[4].set_ylim([0, 100 * 1.05])

    ylim_maxs = {}
    for resolution in convergence_resolutions:
        sim_id = f"lobule_sixth_{resolution}"
        label = convergence_labels[resolution]
        color = convergence_colors[resolution]

        # interpolate time
        xr_cells_raw = xr_cells_dict[sim_id]
        xr_cells = xr_cells_raw.interp(time=times)

        for k_col, sid in enumerate(
            ["rr_(S_ext)", "rr_(P_ext)", "rr_protein", "rr_(T)", "rr_necrosis"]
        ):
            # console.print(xr_cells)
            x = xr_cells.time / 60 / 60  # [s] -> [hr]

            # weighting mean by volume
            cell_volumes = xr_cells.element_volume_point_TPM
            y = xr_cells[sid].weighted(cell_volumes).mean(dim="cell")
            yerr = xr_cells[sid].weighted(cell_volumes).std(dim="cell")

            # y: xr.Dataset = (xr_cells[sid] * cell_volumes).sum(
            #     dim="cell"
            # ) / cell_volumes.sum(dim="cell")

            if sid == "rr_necrosis":
                # [%]
                y = y * 100
                yerr = yerr * 100

            # update max
            if sid not in ylim_maxs:
                ylim_maxs[sid] = 0.0
            if (y + yerr).max() > ylim_maxs[sid]:
                ylim_maxs[sid] = float((y + yerr).max())

            ax = axes[k_col]
            ax.errorbar(
                x=x,
                y=y,
                # yerr=yerr,
                color=color,
                label=label if sid == "rr_necrosis" else "__nolabel__",
                **kwargs,
            )

    for ax in axes:
        ax.set_ylim(bottom=0)
    fig.legend(
        loc="lower center",
        bbox_to_anchor=(0.5, 1.0),
        ncol=n_cols,
        frameon=True,
    )
    plt.tight_layout()
    plt.show()
    fig.savefig(results_dir / "convergence_time.png", bbox_inches="tight")


def plot_convergence_over_position(
    results_dir: Path,
    xr_cells_dict: dict[str, xr.Dataset],
) -> None:
    """Plot convergence over position."""
    console.rule(title="position", style="white")

    kwargs_position = {
        "linestyle": "",
        "marker": "o",
        "markeredgecolor": "black",
        "markeredgewidth": 0.5,
        "markersize": 3,
    }

    n_cols = len(convergence_resolutions)
    fig, axes = plt.subplots(
        ncols=n_cols,
        figsize=(n_cols * 2.7, 2.5),
        dpi=300,
        layout="constrained",
    )

    for ax in axes:
        ax.set_xlabel("Position [-]", **label_kwargs)

    axes[0].set_ylabel("Substrate [mM]", **label_kwargs)
    axes[1].set_ylabel("Product [mM]", **label_kwargs)
    axes[2].set_ylabel("Protein [-]", **label_kwargs)
    axes[3].set_ylabel("Toxic compound [mM]", **label_kwargs)
    axes[4].set_ylabel("Necrosis [%]", **label_kwargs)

    ylim_maxs = {}
    for resolution in convergence_resolutions:
        sim_id = f"lobule_sixth_{resolution}"
        label = convergence_labels[resolution]
        color = convergence_colors[resolution]

        xr_cells_raw = xr_cells_dict[sim_id]

        # interpolate time (only last timepoint)
        xr_cells = xr_cells_raw.interp(time=1.5 * 60 * 60)  # 1.5 hr

        for k_col, sid in enumerate(
            ["rr_(S_ext)", "rr_(P_ext)", "rr_protein", "rr_(T)", "rr_necrosis"]
        ):
            x = xr_cells.rr_position
            y = xr_cells[sid]

            # filter necrosis values != 0.0 or 1.0 (partial necrosis due to point averaging)
            if sid == "rr_necrosis":
                x = x.where((y == 0.0) | (y == 1.0))
                y = y.where((y == 0.0) | (y == 1.0))

                # [%]
                y = y * 100

            # update max
            if sid not in ylim_maxs:
                ylim_maxs[sid] = 0.0
            if y.max() > ylim_maxs[sid]:
                ylim_maxs[sid] = float(y.max())

            axes[k_col].plot(
                x,
                y,
                label=label if sid == "rr_necrosis" else "__nolabel__",
                color=color,
                **kwargs_position,
            )

    for k_col, sid in enumerate(
        ["rr_(S_ext)", "rr_(P_ext)", "rr_protein", "rr_(T)", "rr_necrosis"]
    ):
        axes[k_col].set_ylim([-0.05 * ylim_maxs[sid], 1.05 * ylim_maxs[sid]])

    for ax in axes:
        ax.set_ylim(bottom=0)
        ax.xaxis.set_ticks([0, 0.5, 1], labels=["PP", "", "PV"])

    fig.legend(
        loc="lower center",
        bbox_to_anchor=(0.5, 1.0),
        ncol=n_cols,
        frameon=True,
    )
    plt.tight_layout()
    plt.show()
    fig.savefig(results_dir / "convergence_position.png", bbox_inches="tight")


if __name__ == "__main__":
    """Analysis plots of the convergence simulation."""
    from porous_media.analyses.spt import data_spt_dir

    # XDMF
    xdmf_dir = data_spt_dir / "convergence_sixth" / "xdmf"
    results_dir = data_spt_dir / "convergence_sixth"
    results_dir.mkdir(exist_ok=True, parents=True)

    xdmf_paths = sorted([f for f in xdmf_dir.glob("*.xdmf")])

    # Load xarray datasets
    xr_cells_dict: dict[str, xr.Dataset] = {}
    xr_points_dict: dict[str, xr.Dataset] = {}
    tend: float = np.inf
    for xdmf_path in xdmf_paths:
        xr_cells, xr_points = mesh_datasets_from_xdmf(xdmf_path)

        # get sim id and store data
        sim_id = xdmf_path.stem
        xr_cells_dict[sim_id] = xr_cells
        xr_points_dict[sim_id] = xr_points

        tend_sim: float = xr_cells.time.values[-1]
        if tend_sim < tend:
            tend = tend_sim

    # figure out end time
    time_vec: np.ndarray = np.linspace(start=0, stop=tend, num=101)

    plot_convergence_over_time(
        results_dir=results_dir,
        xr_cells_dict=xr_cells_dict,
        times=time_vec,
    )

    plot_convergence_over_position(
        results_dir=results_dir,
        xr_cells_dict=xr_cells_dict,
    )

# TODO: resolution dependency
