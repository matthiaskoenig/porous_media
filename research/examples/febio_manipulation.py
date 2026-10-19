"""Archived example; configure research paths before running."""

from pathlib import Path

import numpy as np

from porous_media.febio.febio_manipulation import create_feb_files

if __name__ == "__main__":
    from research import DATA_DIR

    feb_template = DATA_DIR / "febio" / "lobule_BCflux.feb.template"

    # -----------------------------
    # Example1: simple rendering
    # -----------------------------
    # console.rule(title="Example1", align="left", style="white")
    print("-" * 80)
    feb_out_dir: Path = DATA_DIR / "febio" / "febs_example1"

    feb_data = {
        "sim1": {"boundary_pressure": 133, "boundary_flux": 1},
        "sim2": {"boundary_pressure": 266, "boundary_flux": 2},
    }
    create_feb_files(
        feb_template=feb_template,
        feb_data=feb_data,
        output_dir=feb_out_dir,
    )

    # --------------------------------
    # Example2: use numpy vectors
    # --------------------------------
    # console.rule(title="Example2", align="left", style="white")
    print("-" * 80)
    boundary_pressure = np.linspace(start=1.37 * 133, stop=5 * 133, num=10)
    boundary_flow = np.linspace(start=0.00008441679, stop=0.00000422083951, num=10)

    feb_data = {}
    counter = 1
    for pressure in boundary_pressure:
        for flow in boundary_flow:
            feb_data[f"sim{counter:>03}"] = {
                "boundary_pressure": pressure,
                "boundary_flow": flow,
            }
            counter += 1

    feb_template = DATA_DIR / "febio" / "lobule_BCflux.feb.template"
    feb_out_dir = DATA_DIR / "febio" / "scanpf"
    create_feb_files(
        feb_template=feb_template,
        feb_data=feb_data,
        output_dir=feb_out_dir,
    )
