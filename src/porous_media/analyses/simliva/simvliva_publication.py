"""Visualization of Simliva results.

Reading results and creating visualizations.
"""
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import meshio

from porous_media import DATA_DIR
from porous_media.console import console
from porous_media.mesh.mesh_tools import MeshTimepoint
from porous_media.visualization.pyvista_visualization import visualize_lobulus_vtk, \
    Scalar, calculate_value_ranges

from porous_media import BASE_DIR
from porous_media.visualization.image_manipulation import merge_images


def figure_temperature_dependency(scalars: List[Scalar]):
    """Temperature dependency."""

    xdmfs: Dict[str, Path] = {
        "sim_T277": DATA_DIR / "simliva" / "005_T_277_15K_P0__0Pa_t_24h" / "results_interpolated_11.xdmf",
        "sim_T310": DATA_DIR / "simliva" / "005_T_277_15K_P0__0Pa_t_24h" / "results_interpolated_11.xdmf",
    }

    # add limits
    for sim_key, xdmf_path in xdmfs.items():
        calculate_value_ranges(xdmf_path=xdmf_path, scalars=scalars)
        console.rule(title="Scalars", align="left", style="white")
        console.print(scalars_iri)



if __name__ == "__main__":

    output_path = BASE_DIR / "results" / "simliva_publication"

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
    exit()

    # Create combined figures for variables and timepoints
    scalars_plot = ["GLC", "O2", "LAC", "ATP", "ADP", "ROS", "necrosis", "ALT", "AST"]
    # create figures for individual panels
    visualize_panels(xdmf_path=xdmf_path, scalars=scalars_iri, output_path=output_path)

    # static images
    # for sim_key, vtk_paths in vtks.items():
    #     row_dir: Path = output_path / sim_key / "rows"
    #     row_dir.mkdir(parents=True, exist_ok=True)
    #
    #     rows: List[Path] = []
    #     for p in vtk_paths:
    #         row: List[Path] = []
    #         for scalar in scalars_plot:
    #             img_path = output_path / sim_key / "panels" / scalar / f"{p.stem}.png"
    #             row.append(img_path)
    #
    #         row_image: Path = output_path / sim_key / "rows" / f"{p.stem}.png"
    #         merge_images(paths=row, direction="horizontal", output_path=row_image)
    #         rows.append(row_image)
    #
    #     console.print(rows)
    #     image: Path = output_path / sim_key / f"{sim_key}.png"
    #     merge_images(paths=rows, direction="vertical", output_path=image)

    # gifs
    for k in range(len(times_wanted)):
        rows: List[Path] = []
        for sim_key, paths in vtks.items():
            row_dir: Path = output_path / sim_key / "rows"
            row_dir.mkdir(parents=True, exist_ok=True)

            vtk_path = paths[k]
            row: List[Path] = []
            for scalar in scalars_plot:
                img_path = output_path / sim_key / "panels" / scalar / f"{p.stem}.png"
                row.append(img_path)

            row_image: Path = output_path / sim_key / "rows" / f"{p.stem}.png"
            print(row)
            merge_images(paths=row, direction="horizontal", output_path=row_image)
            rows.append(row_image)

        console.print(rows)
        image: Path = output_path / sim_key / "gifs" / f"{k:<03}.png"
        merge_images(paths=rows, direction="vertical", output_path=image)
