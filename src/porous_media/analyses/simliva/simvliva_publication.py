"""Visualization of Simliva results."""
from pathlib import Path
from typing import Dict, List

import meshio

from porous_media import DATA_DIR
from porous_media.console import console
from porous_media.visualization.pyvista_visualization import (
    Scalar,
    calculate_value_ranges, visualize_scalars_timecourse
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

    # for num in [10, 200]:
    # for num in [10, ]:
    #     xdmf_path = DATA_DIR / "simliva" / "iri_flux_study_0" / f"results_interpolated_{num}.xdmf"
    #     output_dir = BASE_DIR / "results" / "simliva" / f"gradient_{num}"
    #
    #     # add data limits
    #     calculate_value_ranges(xdmf_path=xdmf_path, scalars=scalars)
    #     console.rule(title="Scalars", align="left", style="white")
    #     console.print(scalars_iri)
    #
    #     # create panels
    #     visualize_scalars_timecourse(
    #         xdmf_path=xdmf_path,
    #         scalars=scalars_iri,
    #         output_dir=output_dir
    #     )

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

    # Create combined figure for timecourse
    create_combined_figure(
        xdmf_path=DATA_DIR / "simliva" / "iri_flux_study_0" / f"results_interpolated_10.xdmf",
        output_dir=BASE_DIR / "results" / "simliva" / f"gradient_10",
        scalars_selection=scalars_selection,
    )

    # FIXME: create gifs


def create_combined_figure(xdmf_path, output_dir, scalars_selection: List[str]):
    """Create combined image of the panels."""
    row_dir: Path = output_dir / "rows"
    row_dir.mkdir(parents=True, exist_ok=True)

    rows: List[Path] = []
    with meshio.xdmf.TimeSeriesReader(xdmf_path) as reader:
        points, cells = reader.read_points_cells()
        for k in range(reader.num_steps):
            row: List[Path] = []
            for scalar_id in scalars_selection:
                img_path = output_dir / "panels" / scalar_id / f"sim_{k:05d}.png"
                row.append(img_path)

            row_image: Path = output_dir / "rows" / f"sim_{k:05d}.png"
            merge_images(paths=row, direction="horizontal", output_path=row_image)
            rows.append(row_image)

    console.print(rows)
    image: Path = output_dir / f"gradient.png"
    merge_images(paths=rows, direction="vertical", output_path=image)

if __name__ == "__main__":
    # Information for plots
    scalars_iri: List[Scalar] = [
        Scalar(sid="rr_necrosis", title="Necrosis (0: alive, 1: death)", colormap="binary"),
        Scalar(sid="rr_(atp)", title="ATP (mM)", colormap="RdBu"),
        Scalar(sid="rr_(glc)", title="Glucose (mM)", colormap="RdBu"),
        Scalar(sid="rr_(lac)", title="Lactate (mM)", colormap="RdBu"),
        Scalar(sid="rr_(o2)", title="Oxygen (mM)", colormap="RdBu"),
        Scalar(sid="rr_(pyr)", title="Pyruvate (mM)", colormap="RdBu"),
        Scalar(sid="rr_(adp)", title="ADP (mM)", colormap="RdBu"),
        Scalar(sid="rr_(nadh)", title="NADH (mM)", colormap="RdBu"),
        Scalar(sid="rr_(nad)", title="NAD (mM)", colormap="RdBu"),
        Scalar(sid="rr_(ros)", title="ROS (mM)", colormap="RdBu"),
        Scalar(sid="rr_(alt)", title="ALT (mM)", colormap="RdBu"),
        Scalar(sid="rr_(ast)", title="AST (mM)", colormap="RdBu"),
        Scalar(sid="rr_Vext", title="Volume extern (l)", colormap="RdBu"),
        Scalar(sid="rr_Vli", title="Volume liver (l)", colormap="RdBu"),
    ]

    # figure_temperature_dependency(scalars=scalars_iri)
    figure_gradient(scalars=scalars_iri)
