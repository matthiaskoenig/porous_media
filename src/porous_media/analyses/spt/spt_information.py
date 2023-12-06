"""Information on the boundary scans"""

from pathlib import Path

import numpy as np
import pandas as pd
from porous_media.console import console
import matplotlib
from matplotlib import pyplot as plt

# substrate boundary flows (2023-12-05)
boundary_flows = [
    -1.0937948e-05,
    -1.2713134e-05,
    -1.448832e-05,
    -1.6263506e-05,
    -1.8038692e-05
]
# substrate boundary flows (2023-12-06)
# boundary_flows = [
#     -1.09379480e-07, -2.93431871e-07, -7.87188447e-07, -2.11178714e-06,
#     -5.66528247e-06, -1.51982295e-05, -4.07722263e-05, -1.09379480e-04
# ]

# zonation patterns
pattern_idx2name = {
    0: "constant",
    1: "linear_increase",
    2: "linear decrease",
    3: "sharp_pericentral",
    4: "sharp_periportal",
}
pattern_name2idx = {v: k for k, v in pattern_idx2name.items()}

pattern_order = [
    "constant",
    "linear_increase",
    "sharp_pericentral",
    "linear decrease",
    "sharp_periportal",
]

# colors
pattern_colors = {
    0: "tab:blue",
    1: "tab:orange",
    2: "tab:green",
    3: "tab:red",
    4: "tab:purple",
}
pattern_colormaps = {
    0: "Greys",
    1: "Blues",
    2: "Greens",
    3: "Reds",
    4: "Purples",
}


def simulation_conditions_df() -> pd.DataFrame:
    """Create simulation condition DataFrame."""

    feb_data = {}
    counter = 1
    for pattern_key, pattern_name in pattern_idx2name.items():
        kmax = len(boundary_flows)
        for k_flow, flow in enumerate(boundary_flows):

            # calculate colors
            cmap_key = pattern_colormaps[pattern_key]
            cmap = matplotlib.cm.get_cmap(cmap_key)
            # color

            # dose_relative = np.log(dose) / np.log(dose_max)
            color_rgba = cmap((k_flow + 1) / kmax)
            color_hex = matplotlib.colors.to_hex(color_rgba, keep_alpha=True)

            sim_id = f"sim{counter:>03}"

            feb_data[sim_id] = {
                "pattern_key": pattern_key,
                "pattern_name": pattern_name,
                "boundary_flow_key": k_flow,
                "boundary_flow": flow,
                "color": color_hex,
            }
            counter += 1

    df = pd.DataFrame(feb_data).T
    df = df.astype({"pattern_key": int, "boundary_flow_key": int, "color": str})
    console.print(df)
    return df

plot_kwargs = {
    "marker": "o",
    "markersize": 10,
    "markeredgecolor": "black",
}

def plot_boundary_flux(df: pd.DataFrame):
    """Analysis of boundary scan"""
    y = np.abs(df.boundary_flow.values[0:len(boundary_flows)])

    f, ax = plt.subplots(nrows=1, ncols=1, figsize=(5, 5), dpi=300, layout="constrained")

    ax.plot(
        y,
        label="simulation",
        **plot_kwargs,
    )

    ynew = np.linspace(start=y[-1]*0.01, stop=y[-1], num=8)
    console.print(ynew)
    ax.plot(
        ynew,
        label="ynew",
        **plot_kwargs,
    )

    ax.set_ylim(bottom=0)
    ax.set_xlabel("simulation index", fontsize=15, fontweight="bold")
    ax.set_ylabel("absolute boundary flow", fontsize=15, fontweight="bold")
    ax.legend()
    plt.show()


def plot_colors(df: pd.DataFrame):
    """Plot colors"""
    f, ax = plt.subplots(nrows=1, ncols=1, figsize=(5, 5), dpi=300, layout="constrained")

    for sim_id, row in df.iterrows():
        ax.plot(
            row.pattern_key,
            row.boundary_flow_key,
            linestyle=None,
            color=row.color,
            **plot_kwargs,
        )

    ax.set_ylim(bottom=0)
    ax.set_xlabel("zonation pattern", fontsize=15, fontweight="bold")
    ax.set_ylabel("boundary flow", fontsize=15, fontweight="bold")
    plt.show()


if __name__ == "__main__":
    df = simulation_conditions_df()
    # xdmf_dir = Path("/home/mkoenig/git/porous_media/data/spt/2023-12-05")
    xdmf_dir = Path("/home/mkoenig/git/porous_media/data/spt/2023-12-06")
    df.to_excel(xdmf_dir / "information.xlsx", index=True)

    plot_boundary_flux(df)
    plot_colors(df)
