"""Substrate scan simulation."""

from pathlib import Path
from typing import Dict, List

import numpy as np

from porous_media import BASE_DIR, DATA_DIR
from porous_media.console import console
from porous_media.febio.results_processing import interpolate_xdmf, vtks_to_xdmf, \
    xdmfs_from_febio, xdmf_information, XDMFInformation
from porous_media.visualization.image_manipulation import merge_images
from porous_media.visualization.pyvista_visualization import (
    Scalar,
    calculate_value_ranges,
    create_combined_images,
    visualize_scalars_timecourse,
)
from porous_media.visualization.video import create_video
from porous_media.log import get_logger

logger = get_logger(__name__)


def visualize_substrate_videos(xdmf_paths: List[Path], scalars: List[Scalar], results_dir: Path, create_panels: bool = True) -> None:
    """Videos for substrate dependency"""

    tend = 28800
    if create_panels:
        for xdmf_path in xdmf_paths:

            for num in [10, 250]:
                output_dir = results_dir / f"{xdmf_path.stem}_{num}"

                # interpolate
                interpolate_xdmf(
                    xdmf_in=xdmf_path,
                    xdmf_out=xdmf_path.parent / f"{xdmf_path.stem}_interpolated_{num}.xdmf",
                    times_interpolate=np.linspace(0, tend, num=num),
                )

                # add data limits
                # FIXME: this should be independent of scalars
                calculate_value_ranges(xdmf_path=xdmf_path, scalars=scalars)
                console.rule(title="Scalars", align="left", style="white")
                console.print(scalars_spt)

                # create panels
                visualize_scalars_timecourse(
                    xdmf_path=xdmf_path, scalars=scalars_spt, output_dir=output_dir
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
    xdmf_dir = Path("/home/mkoenig/git/porous_media/data/spt_substrate_scan")
    xdmf_paths: List[Path] = xdmfs_from_febio(
        # febio_dir=Path("/home/mkoenig/data/qualiperf/P7-Perf/spt_results/simulation_fixedcelltype"),
        febio_dir=Path("/home/mkoenig/git/porous_media/data/spt/simulation_fixedcelltype"),
        xdmf_dir=xdmf_dir,
        overwrite=False,
    )

    xdmf_paths = [xdmf_dir / f"sim_{k}.xdmf" for k in range(21, 26)]
    console.print(xdmf_paths)
    results_dir: Path = BASE_DIR / "results" / "spt_substrate_scan"

    info: XDMFInformation = xdmf_information(xdmf_path=xdmf_paths[0])
    console.print(info)

    # from porous_media.analyses.spt import scalars_spt
    # visualize_substrate_videos(
    #     xdmf_paths=xdmf_paths,
    #     scalars=scalars_spt,
    #     results_dir=results_dir,
    # )
