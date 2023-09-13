"""Visualization of Simliva results."""
from pathlib import Path
from typing import Dict, List

from porous_media import DATA_DIR
from porous_media.console import console
from porous_media.visualization.pyvista_visualization import (
    Scalar,
    calculate_value_ranges,
)

from porous_media import BASE_DIR
from porous_media.visualization.image_manipulation import merge_images


def figure_temperature_dependency(scalars: List[Scalar]):
    """Temperature dependency."""
    output_path = BASE_DIR / "results" / "simliva" / "temperature_dependency"
    output_path.mkdir(exist_ok=True, parents=True)

    xdmfs: Dict[str, Path] = {
        "sim_T277": DATA_DIR / "simliva" / "005_T_277_15K_P0__0Pa_t_24h" / "results_interpolated_11.xdmf",
        "sim_T310": DATA_DIR / "simliva" / "005_T_277_15K_P0__0Pa_t_24h" / "results_interpolated_11.xdmf",
    }

    # add data limits
    for sim_key, xdmf_path in xdmfs.items():
        calculate_value_ranges(xdmf_path=xdmf_path, scalars=scalars)
        console.rule(title="Scalars", align="left", style="white")
        console.print(scalars_iri)

    # create panels for all timepoints and scalars
    visualize_panels(xdmf_path=xdmf_path, scalars=scalars_iri, output_path=output_path)


    # Create combined figure for all timecourses
    scalars_plot = ["GLC", "O2", "LAC", "ATP", "ADP", "ROS", "necrosis", "ALT", "AST"]
    for sim_key, vtk_paths in vtks.items():
        row_dir: Path = output_path / sim_key / "rows"
        row_dir.mkdir(parents=True, exist_ok=True)

        rows: List[Path] = []
        for p in vtk_paths:
            row: List[Path] = []
            for scalar in scalars_plot:
                img_path = output_path / sim_key / "panels" / scalar / f"{p.stem}.png"
                row.append(img_path)

            row_image: Path = output_path / sim_key / "rows" / f"{p.stem}.png"
            merge_images(paths=row, direction="horizontal", output_path=row_image)
            rows.append(row_image)

        console.print(rows)
        image: Path = output_path / sim_key / f"{sim_key}.png"
        merge_images(paths=rows, direction="vertical", output_path=image)


def figure_gradient(scalars: List[Scalar]):
    """Gradient dependency."""
    output_path = BASE_DIR / "results" / "simliva" / "gradient"
    output_path.mkdir(exist_ok=True, parents=True)
    xdmf_path = DATA_DIR / "simliva" / "iri_flux_study_0" / "results_interpolated.xdmf"

    # add data limits
    calculate_value_ranges(xdmf_path=xdmf_path, scalars=scalars)
    console.rule(title="Scalars", align="left", style="white")
    console.print(scalars_iri)

    # create panels for all timepoints and scalars
    visualize_panels(xdmf_path=xdmf_path, scalars=scalars_iri, output_path=output_path)



    # Create combined figure for all timecourses
    scalars_plot = ["GLC", "O2", "LAC", "ATP", "ADP", "ROS", "necrosis", "ALT", "AST"]
    for sim_key, vtk_paths in vtks.items():
        row_dir: Path = output_path / sim_key / "rows"
        row_dir.mkdir(parents=True, exist_ok=True)

        rows: List[Path] = []
        for p in vtk_paths:
            row: List[Path] = []
            for scalar in scalars_plot:
                img_path = output_path / sim_key / "panels" / scalar / f"{p.stem}.png"
                row.append(img_path)

            row_image: Path = output_path / sim_key / "rows" / f"{p.stem}.png"
            merge_images(paths=row, direction="horizontal", output_path=row_image)
            rows.append(row_image)

        console.print(rows)
        image: Path = output_path / sim_key / f"{sim_key}.png"
        merge_images(paths=rows, direction="vertical", output_path=image)


if __name__ == "__main__":
    # Information for plots
    scalars_iri: List[Scalar] = [
        Scalar(sid="necrosis", title="Necrosis (0: alive, 1: death)", cmap="binary"),
        Scalar(sid="ATP", title="ATP (mM)", cmap="RdBu"),
        Scalar(sid="GLC", title="Glucose (mM)", cmap="RdBu"),
        Scalar(sid="LAC", title="Lactate (mM)", cmap="RdBu"),
        Scalar(sid="O2", title="Oxygen (mM)", cmap="RdBu"),
        Scalar(sid="PYR", title="Pyruvate (mM)", cmap="RdBu"),
        Scalar(sid="ADP", title="ADP (mM)", cmap="RdBu"),
        Scalar(sid="NADH", title="NADH (mM)", cmap="RdBu"),
        Scalar(sid="NAD", title="NAD (mM)", cmap="RdBu"),
        Scalar(sid="ROS", title="ROS (mM)", cmap="RdBu"),
        Scalar(sid="ALT", title="ALT (mM)", cmap="RdBu"),
        Scalar(sid="AST", title="AST (mM)", cmap="RdBu"),
    ]

    figure_temperature_dependency(scalars=scalars_iri)
