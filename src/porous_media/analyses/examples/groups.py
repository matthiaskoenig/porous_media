"""Visualization of example results."""
from pathlib import Path
from typing import Dict, List

import numpy as np

from porous_media import BASE_DIR, DATA_DIR
from porous_media.console import console
from porous_media.data.xdmf_tools import interpolate_xdmf, vtks_to_xdmf
from porous_media.visualization.image_manipulation import merge_images
from porous_media.visualization.pyvista_visualization import (
    DataLayer,
    create_combined_images,
    visualize_datalayers_timecourse,
)
from porous_media.visualization.video import create_video


def visualize_group(scalars: List[DataLayer], create_panels: bool = True) -> None:
    """Gradient dependency."""
    sim_dir = DATA_DIR / "examples" / "Group_PF_ideal_TPM"
    xdmf_path = sim_dir / "results.xdmf"

    num_steps = 100
    xdmf_interpolated_path = sim_dir / f"results_interpolated_{num_steps}.xdmf"
    output_dir = BASE_DIR / "results" / "groups"

    if create_panels:
        # process VTKs
        vtks_to_xdmf(vtk_dir=sim_dir / "vtk", xdmf_path=xdmf_path)

        # interpolate
        interpolate_xdmf(
            xdmf_in=xdmf_path,
            xdmf_out=xdmf_interpolated_path,
            times_interpolate=np.linspace(0, 1, num=num_steps),
        )

        # add data limits

        console.rule(title="Scalars", align="left", style="white")
        console.print(scalars)

        # create panels
        visualize_datalayers_timecourse(
            xdmf_path=xdmf_interpolated_path, data_layers=scalars, output_dir=output_dir
        )

    # subset of scalars to visualize
    scalars_selection: List[str] = [
        "stress",
        "fluid_flux_TPM",
        "seepage_velocity_TPM",
        "solid_stress_TPM",
        "Lagrange_strain",
        "pressure",
    ]

    # Create combined image
    rows: List[Path] = create_combined_images(
        num_steps=num_steps,
        output_dir=output_dir,
        selection=scalars_selection,
        direction="horizontal",
    )
    rows_subset = [rows[k] for k in range(len(rows)) if (k + 1) % 10 == 0]
    merge_images(
        paths=rows_subset, direction="vertical", output_path=output_dir / "groups.png"
    )

    # Create video
    create_combined_images(
        num_steps=num_steps,
        output_dir=output_dir,
        selection=scalars_selection,
        direction="custom",
        ncols=3,
        nrows=2,
    )
    create_video(
        image_pattern=str(output_dir / "custom" / "sim_%05d.png"),
        video_path=output_dir / "groups.mp4",
    )


if __name__ == "__main__":
    scalars: List[DataLayer] = [
        DataLayer(sid="stress", title="Stress", colormap="RdBu"),
        DataLayer(sid="fluid_flux_TPM", title="Fluid flux TPM", colormap="RdBu"),
        DataLayer(
            sid="seepage_velocity_TPM", title="Seepage velocity TPM", colormap="RdBu"
        ),
        DataLayer(sid="solid_stress_TPM", title="Solid stress", colormap="RdBu"),
        DataLayer(sid="Lagrange_strain", title="Lagrange strain", colormap="RdBu"),
        DataLayer(sid="pressure", title="Pressure", colormap="RdBu"),
    ]
    visualize_group(scalars=scalars, create_panels=False)
