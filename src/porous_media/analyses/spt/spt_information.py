"""Information on the boundary scans"""

from pathlib import Path

import numpy as np
import pandas as pd
from porous_media.console import console
from matplotlib import pyplot as plt

pattern_names = {
    0: "constant",
    1: "linear_increase",
    2: "linear_decrease",
    3: "sharp_pericentral",
    4: "sharp_periportal",
}

pattern_colors = {
    0: "tab:blue",
    1: "tab:orange",
    2: "tab:green",
    3: "sharp_pericentral",
    4: "sharp_periportal",
}

colors = [

    "tab:orange",
    "tab:green",
    "tab:red",
    "tab:purple",
]

def simulation_conditions_df() -> pd.DataFrame:
    """Create simulation condition DataFrame."""



    zonation_pattern = [0, 1, 2, 3, 4]
    boundary_flow = [-1.0937948e-05, -1.2713134e-05, -1.448832e-05, -1.6263506e-05, -1.8038692e-05]

    feb_data = {}
    counter = 1
    for pattern in zonation_pattern:
        for flow in boundary_flow:
            feb_data[f"sim{counter:>03}"] = {
                "pattern": pattern,
                "pattern_name": pattern_names[pattern],
                "boundary_flow": flow,
            }
            counter += 1

    df = pd.DataFrame(feb_data).T
    df = df.astype({"pattern": int})
    console.print(df)
    return df


def plot_boundary_flux(df: pd.DataFrame):
    y = df.boundary_flow.values

    f, ax = plt.subplots(nrows=1, ncols=1)
    ax.plot(
        np.abs(y[0:5]),
        marker="o",
        label="simulation"
    )

    ynew = np.logspace(-2, 1, num=8)*y[5]
    console.print(ynew)
    ax.plot(
        np.abs(ynew),
        marker="o",
        label="new scan"
    )

    ax.set_ylim(bottom=0)
    ax.set_xlabel("zonation pattern")
    ax.set_ylabel("absolute boundary flow")

    plt.show()


if __name__ == "__main__":
    df = simulation_conditions_df()
    xdmf_dir = Path("/home/mkoenig/git/porous_media/data/spt/2023-12-05")
    df.to_excel(xdmf_dir / "simulations.xlsx", index=True)

    plot_boundary_flux(df)
