"""Substrate scan simulation."""

from pathlib import Path
from typing import Dict, List

import numpy as np
from porous_media.console import console
from porous_media.febio.results_processing import xdmfs_from_febio
from porous_media.visualization.image_manipulation import merge_images
from porous_media.visualization.pyvista_visualization import (
    Scalar,
    calculate_value_ranges,
    create_combined_images,
    visualize_scalars_timecourse,
)
from porous_media.visualization.video import create_video
from porous_media.analyses.spt import scalars_spt
from porous_media.log import get_logger

logger = get_logger(__name__)

if __name__ == "__main__":

    xdmf_paths: List[Path] = xdmfs_from_febio(
        febio_dir=Path("/home/mkoenig/data/qualiperf/P7-Perf/spt_results/simulation_zonation_hard"),
        xdmf_dir=Path("/home/mkoenig/git/porous_media/data/spt_zonation_patterns"),
        overwrite=False,
    )

    # visualize_gradient(scalars=scalars_iri, create_panels=False)
