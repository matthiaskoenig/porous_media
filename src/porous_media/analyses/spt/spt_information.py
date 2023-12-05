"""Information on the boundary scans"""

from pathlib import Path

import numpy as np
import pandas as pd
from porous_media.console import console
import matplotlib
from matplotlib import pyplot as plt

# substrate boundary flows
boundary_flows = [
    -1.0937948e-05,
    -1.2713134e-05,
    -1.448832e-05,
    -1.6263506e-05,
    -1.8038692e-05
]

# zonation patterns
pattern_names = {
    0: "constant",
    1: "linear_increase",
    2: "linear_decrease",
    3: "sharp_pericentral",
    4: "sharp_periportal",
}

# colors
pattern_colors = {
    0: "tab:blue",
    1: "tab:orange",
    2: "tab:green",
    3: "tab:red",
    4: "tab:purple",
}
pattern_colormaps = {
    0: "Blues",
    1: "Oranges",
    2: "Greens",
    3: "Reds",
    4: "Purples",
}




def simulation_conditions_df() -> pd.DataFrame:
    """Create simulation condition DataFrame."""

    zonation_pattern = [0, 1, 2, 3, 4]

    feb_data = {}
    counter = 1
    for pattern_key in zonation_pattern:
        kmax = len(boundary_flows)
        for k_flow, flow in enumerate(boundary_flows):

            # calculate colors
            cmap_key = pattern_colormaps[pattern_key]
            cmap = matplotlib.cm.get_cmap(cmap_key)
            # color

            # dose_relative = np.log(dose) / np.log(dose_max)
            color_rgba = cmap((k_flow + 1) / kmax)
            color_hex = matplotlib.colors.to_hex(color_rgba, keep_alpha=True)

            # sim_id = f"sim{counter:>03}"
            sim_id = f"sim00{counter}"  # FIXME: bug by Steffen in ids

            feb_data[sim_id] = {
                "pattern_key": pattern_key,
                "pattern_name": pattern_names[pattern_key],
                "boundary_flow_key": k_flow,
                "boundary_flow": flow,
                "color": color_hex,
            }
            counter += 1

    df = pd.DataFrame(feb_data).T
    df = df.astype({"pattern_key": int, "boundary_flow_key": int, "color": str})
    console.print(df)
    return df


def plot_boundary_flux(df: pd.DataFrame):
    """Analysis of boundary scan"""
    y = df.boundary_flow.values

    f, ax = plt.subplots(nrows=1, ncols=1)
    ax.plot(
        np.abs(y[0:5]),
        marker="o",
        label="simulation"
    )

    ynew = np.logspace(-2, 1, num=8)*y[5]
    # console.print(ynew)
    ax.plot(
        np.abs(ynew),
        marker="o",
        label="new scan"
    )

    ax.set_ylim(bottom=0)
    ax.set_xlabel("simulation index")
    ax.set_ylabel("absolute boundary flow")

    plt.show()


if __name__ == "__main__":
    df = simulation_conditions_df()
    xdmf_dir = Path("/home/mkoenig/git/porous_media/data/spt/2023-12-05")
    df.to_excel(xdmf_dir / "information.xlsx", index=True)

    plot_boundary_flux(df)
