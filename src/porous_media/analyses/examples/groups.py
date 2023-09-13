"""Visualization of example results."""
from pathlib import Path
from typing import Dict, List

import numpy as np

from porous_media import BASE_DIR, DATA_DIR
from porous_media.console import console
from porous_media.febio.results_processing import interpolate_xdmf, vtks_to_xdmf
from porous_media.visualization.image_manipulation import merge_images
from porous_media.visualization.pyvista_visualization import (
    Scalar,
    calculate_value_ranges,
    create_combined_images,
    visualize_scalars_timecourse,
)
from porous_media.visualization.video import create_video


def visualize_group(scalars: List[Scalar], create_panels: bool = True) -> None:
    """Gradient dependency."""
    sim_dir = DATA_DIR / "examples" / "Group_PF_ideal_TPM"
    xdmf_path = sim_dir / f"results.xdmf"
    xdmf_interpolated_path = sim_dir / f"results_interpolated_100.xdmf"
    output_dir = BASE_DIR / "results" / "groups"

    if create_panels:
        # process VTKs
        vtks_to_xdmf(vtk_dir=sim_dir / "vtk", xdmf_path=xdmf_path)

        # interpolate
        interpolate_xdmf(
            xdmf_in=xdmf_path,
            xdmf_out=xdmf_interpolated_path,
            times_interpolate=np.linspace(0, 1, num=100),
        )

        # add data limits
        calculate_value_ranges(xdmf_path=xdmf_interpolated_path, scalars=scalars)
        console.rule(title="Scalars", align="left", style="white")
        console.print(scalars)

        # create panels
        visualize_scalars_timecourse(
            xdmf_path=xdmf_interpolated_path, scalars=scalars, output_dir=output_dir
        )

    # subset of scalars to visualize
    scalars_selection: List[str] = [
        'stress',
        'fluid_flux_TPM',
        'seepage_velocity_TPM',
        'solid_stress_TPM',
        'Lagrange_strain',
        'pressure',
    ]

    # Create combined image
    rows: List[Path] = create_combined_images(
        xdmf_path=xdmf_interpolated_path,
        output_dir=output_dir,
        scalars_selection=scalars_selection,
        direction="horizontal",
    )
    rows_subset = [rows[k] for k in range(len(rows)) if (k+1) % 10 == 0]
    merge_images(
        paths=rows_subset, direction="vertical", output_path=output_dir / "groups.png"
    )

    # Create video
    create_combined_images(
        xdmf_path=xdmf_interpolated_path,
        output_dir=output_dir,
        scalars_selection=scalars_selection,
        direction="custom",
        ncols=3,
        nrows=2,
    )
    create_video(
        image_pattern=str(output_dir / "custom" / "sim_%05d.png"),
        video_path=output_dir / "groups.mp4",
    )


if __name__ == "__main__":
    scalars: List[Scalar] = [
        Scalar(sid="stress", title="Stress", colormap="RdBu"),
        Scalar(sid="fluid_flux_TPM", title="Fluid flux TPM", colormap="RdBu"),
        Scalar(sid="seepage_velocity_TPM", title="Seepage velocity TPM", colormap="RdBu"),
        Scalar(sid="solid_stress_TPM", title="Solid stress", colormap="RdBu"),
        Scalar(sid="Lagrange_strain", title="Lagrange strain", colormap="RdBu"),
        Scalar(sid="pressure", title="Pressure", colormap="RdBu"),
    ]
    visualize_group(scalars=scalars, create_panels=False)
