"""Visualization of Simliva results."""
from pathlib import Path
from typing import Dict, List

import numpy as np

from porous_media import BASE_DIR, DATA_DIR
from porous_media.console import console
from porous_media.febio.xdmf_tools import interpolate_xdmf
from porous_media.visualization.image_manipulation import merge_images
from porous_media.visualization.pyvista_visualization import (
    DataLayer,
    calculate_value_ranges,
    create_combined_images,
    visualize_datalayers_timecourse,
)
from porous_media.visualization.video import create_video


# def visualize_temperature_dependency(scalars: List[Scalar], create_panels: bool = True):
#     """Temperature dependency."""
#     scalars_plot = ["GLC", "O2", "LAC", "ATP", "ADP", "ROS", "necrosis", "ALT", "AST"]
#     output_path = BASE_DIR / "results" / "simliva" / "temperature_dependency"
#     output_path.mkdir(exist_ok=True, parents=True)
#
#     xdmfs: Dict[str, Path] = {
#         "sim_T277": DATA_DIR / "simliva" / "005_T_277_15K_P0__0Pa_t_24h" / "results_interpolated_11.xdmf",
#         "sim_T310": DATA_DIR / "simliva" / "005_T_277_15K_P0__0Pa_t_24h" / "results_interpolated_11.xdmf",
#     }


def visualize_gradient(scalars: List[DataLayer], create_panels: bool = True) -> None:
    """Gradient dependency."""
    sim_dir = DATA_DIR / "simliva" / "iri_flux_study_0"

    if create_panels:
        for num in [10, 250]:
            xdmf_path = sim_dir / f"results_interpolated_{num}.xdmf"
            output_dir = BASE_DIR / "results" / "simliva" / f"gradient_{num}"

            # interpolate
            interpolate_xdmf(
                xdmf_in=sim_dir / "results.xdmf",
                xdmf_out=sim_dir / f"results_interpolated_{num}.xdmf",
                times_interpolate=np.linspace(0, 10000, num=num),
            )

            # add data limits
            calculate_value_ranges(xdmf_path=xdmf_path, scalars=scalars)
            console.rule(title="Scalars", align="left", style="white")
            console.print(scalars_iri)

            # create panels
            visualize_datalayers_timecourse(
                xdmf_path=xdmf_path, data_layers=scalars_iri, output_dir=output_dir
            )

    # subset of scalars to visualize
    scalars_selection: List[str] = [
        "rr_(glc)",
        "rr_(o2)",
        "rr_(lac)",
        "rr_(atp)",
        "rr_(adp)",
        "rr_(ros)",
        "rr_necrosis",
        "rr_(alt)",
        "rr_(ast)",
    ]

    # Create combined images
    for num in [10, 250]:
        output_dir = BASE_DIR / "results" / "simliva" / f"gradient_{num}"
        rows: List[Path] = create_combined_images(
            xdmf_path=DATA_DIR
            / "simliva"
            / "iri_flux_study_0"
            / f"results_interpolated_{num}.xdmf",
            direction="horizontal",
            output_dir=output_dir,
            scalars_selection=scalars_selection,
        )

    # Create combined figure for timecourse
    output_dir = BASE_DIR / "results" / "simliva" / "gradient_10"
    merge_images(
        paths=rows, direction="vertical", output_path=output_dir / "gradient.png"
    )

    # Create video
    output_dir = BASE_DIR / "results" / "simliva" / "gradient_250"
    create_video(
        image_pattern=str(output_dir / "squares" / "sim_%05d.png"),
        video_path=output_dir / "gradient_square.mp4",
    )


if __name__ == "__main__":
    # FIXME: update code for the old simliva example
    # scalars_iri_old: List[Scalar] = [
    #     Scalar(sid="necrosis", title="Necrosis (0: alive, 1: death)",
    #            colormap="binary"),
    #     Scalar(sid="ATP", title="ATP (mM)", colormap="RdBu"),
    #     Scalar(sid="GLC", title="Glucose (mM)", colormap="RdBu"),
    #     Scalar(sid="LAC", title="Lactate (mM)", colormap="RdBu"),
    #     Scalar(sid="O2", title="Oxygen (mM)", colormap="RdBu"),
    #     Scalar(sid="PYR", title="Pyruvate (mM)", colormap="RdBu"),
    #     Scalar(sid="ADP", title="ADP (mM)", colormap="RdBu"),
    #     Scalar(sid="NADH", title="NADH (mM)", colormap="RdBu"),
    #     Scalar(sid="NAD", title="NAD (mM)", colormap="RdBu"),
    #     Scalar(sid="ROS", title="ROS (mM)", colormap="RdBu"),
    #     Scalar(sid="ALT", title="ALT (mM)", colormap="RdBu"),
    #     Scalar(sid="AST", title="AST (mM)", colormap="RdBu"),
    # ]
    # visualize_temperature_dependency(scalars=scalars_iri_old, create_panels=True)

    scalars_iri: List[DataLayer] = [
        DataLayer(
            sid="rr_necrosis",
            title="Necrosis (0: alive, 1: death)",
            colormap="binary",
        ),
        DataLayer(sid="rr_(atp)", title="ATP (mM)", colormap="RdBu"),
        DataLayer(sid="rr_(glc)", title="Glucose (mM)", colormap="RdBu"),
        DataLayer(sid="rr_(lac)", title="Lactate (mM)", colormap="RdBu"),
        DataLayer(sid="rr_(o2)", title="Oxygen (mM)", colormap="RdBu"),
        DataLayer(sid="rr_(pyr)", title="Pyruvate (mM)", colormap="RdBu"),
        DataLayer(sid="rr_(adp)", title="ADP (mM)", colormap="RdBu"),
        DataLayer(sid="rr_(nadh)", title="NADH (mM)", colormap="RdBu"),
        DataLayer(sid="rr_(nad)", title="NAD (mM)", colormap="RdBu"),
        DataLayer(sid="rr_(ros)", title="ROS (mM)", colormap="RdBu"),
        DataLayer(sid="rr_(alt)", title="ALT (mM)", colormap="RdBu"),
        DataLayer(sid="rr_(ast)", title="AST (mM)", colormap="RdBu"),
        DataLayer(sid="rr_Vext", title="Volume extern (l)", colormap="RdBu"),
        DataLayer(sid="rr_Vli", title="Volume liver (l)", colormap="RdBu"),
    ]
    visualize_gradient(scalars=scalars_iri, create_panels=False)
