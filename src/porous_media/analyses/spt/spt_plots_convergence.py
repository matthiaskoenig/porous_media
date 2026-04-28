"""Analysis of convergence."""

from pathlib import Path

import matplotlib
import numpy as np
import xarray as xr
from matplotlib import pyplot as plt

from porous_media.analyses.liver_variables import calculate_necrosis_fraction
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


def plot_spt_over_time(
    results_dir: Path,
    xr_cells_dict: dict[str, xr.Dataset],
    times: np.ndarray,
) -> None:
    """Plot SPT over time."""
    console.rule(title="SPT timecourse", style="white")

    # DataFrame information
    n_cols = 5

    fig, axes = plt.subplots(
        nrows=1,
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
    n_resolutions = len(convergence_resolutions)
    for kres, resolution in enumerate(convergence_resolutions):
        sim_id = f"lobule_sixth_{resolution}"
        label = f"{resolution[4:]}"

        # calculate colors
        cmap = matplotlib.colormaps.get_cmap("viridis")
        color_rgba = cmap(kres / (n_resolutions - 1))
        color = matplotlib.colors.to_hex(color_rgba, keep_alpha=True)

        # interpolate time
        xr_cells_raw = xr_cells_dict[sim_id]
        xr_cells = xr_cells_raw.interp(time=times)

        kwargs = {
            "linestyle": "-",
            "marker": "None",
            "color": color,
            "markeredgecolor": "black",
            "markeredgewidth": 0.5,
            "markersize": 5,
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
                ylim_maxs[sid] = float((y + yerr).max())

            ax = axes[k_col]
            ax.errorbar(
                x=x,
                y=y,
                # yerr=yerr,  FIXME
                # label=label,
                **kwargs,
            )
            # ax.legend()

        necrosis_fraction = calculate_necrosis_fraction(xr_cells=xr_cells)
        axes[4].plot(
            # convert to hr and percent
            necrosis_fraction.time / 60 / 60,  # [s] -> [hr]
            necrosis_fraction * 100,
            label=label,
            **kwargs,
        )

    # for kax, sid in enumerate(["rr_(S_ext)", "rr_(P_ext)", "rr_protein", "rr_(T)"]):
    #     axes[kax].set_ylim([-0.05 * ylim_maxs[sid], 1.05 * ylim_maxs[sid]])

    # for ax in axes[:]:
    #     ax.set_xticks([0, 10, 20])
    #     ax.set_xticklabels([])
    # for ax in axes:
    #     ax.xaxis.set_ticks([0, 10, 20], labels=["0", "10", "20"])

    for ax in axes:
        ax.set_ylim(bottom=0)
    # axes[4].legend()

    fig.legend(
        loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=n_resolutions, frameon=True
    )
    plt.tight_layout()

    plt.show()
    fig.savefig(results_dir / "convergence.png", bbox_inches="tight")


if __name__ == "__main__":
    """Analysis plots of the convergence simulation."""
    from porous_media.analyses.spt import data_spt_dir, results_spt_dir

    # XDMF
    xdmf_dir = Path(data_spt_dir / "convergence_sixth" / "xdmf")
    results_dir = results_spt_dir / "convergence_sixth"
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

    plot_spt_over_time(
        results_dir=results_dir,
        xr_cells_dict=xr_cells_dict,
        times=time_vec,
    )
