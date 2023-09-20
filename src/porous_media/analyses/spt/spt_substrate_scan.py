"""Substrate scan simulation."""

from pathlib import Path
from typing import Dict, List, Set

import numpy as np

from porous_media import BASE_DIR, DATA_DIR
from porous_media.console import console
from porous_media.data.xdmf_tools import (
    DataLimits,
    XDMFInformation,
    interpolate_xdmf,
    xdmfs_from_febio,
)
from porous_media.log import get_logger
from porous_media.visualization.image_manipulation import merge_images
from porous_media.visualization.pyvista_visualization import (
    DataLayer,
    create_combined_images,
    visualize_datalayers_timecourse,
)
from porous_media.visualization.video import create_gif_from_video, create_video


logger = get_logger(__name__)


def visualize_scan(
    xdmf_paths: List[Path],
    data_layers: List[DataLayer],
    results_dir: Path,
    create_panels: bool = True,
) -> None:
    """Videos for substrate dependency"""

    # subset of scalars to visualize
    selection: List[str] = [
        "rr_(S)",
        "rr_(P)",
        "rr_(T)",
        "rr_protein",
        "rr_necrosis",
        "pressure",
    ]

    # Calculate tend time from all simulations
    tends: np.ndarray = np.zeros(shape=(len(xdmf_paths),))
    for k, xdmf_path in enumerate(xdmf_paths):
        xdmf_info = XDMFInformation.from_path(xdmf_path)
        tends[k] = xdmf_info.tend
    tend = tends.min()

    # ordered data layers selected
    data_layers_dict = {dl.sid: dl for dl in data_layers}
    data_layers_selected = [data_layers_dict[sid] for sid in selection]

    if create_panels:
        for num in [10, 100]:
            output_dir = results_dir / f"{num}_{tend}"
            output_dir.mkdir(exist_ok=True, parents=True)

            # interpolate and global limits
            all_limits: List[DataLimits] = []
            for xdmf_path in xdmf_paths:
                # interpolate & calculate limits
                interpolate_xdmf(
                    xdmf_in=xdmf_path,
                    xdmf_out=output_dir / f"{xdmf_path.stem}_interpolated.xdmf",
                    times_interpolate=np.linspace(0, tend, num=num),
                    overwrite=False,
                )

                limits = DataLimits.from_xdmf(xdmf_path=xdmf_path, overwrite=False)
                all_limits.append(limits)

            data_limits = DataLimits.merge_limits(all_limits)
            for data_layer in data_layers_selected:
                data_layer.update_color_limits(data_limits=data_limits)

            # create all panels for interpolation
            for xdmf_path in xdmf_paths:
                visualize_datalayers_timecourse(
                    xdmf_path=output_dir / f"{xdmf_path.stem}_interpolated.xdmf",
                    data_layers=data_layers_selected,
                    output_dir=output_dir / f"{xdmf_path.stem}",
                )

    # Create combined images for all simulations
    for num in [10, 100]:
        for xdmf_path in xdmf_paths:
            output_dir = results_dir / f"{num}_{tend}" / f"{xdmf_path.stem}"
            rows: List[Path] = create_combined_images(
                num_steps=num,
                direction="horizontal",
                output_dir=output_dir,
                selection=selection,
            )

            # Create combined figure for timecourse
            if num == 10:
                merge_images(
                    paths=rows,
                    direction="vertical",
                    output_path=results_dir / f"{xdmf_path.stem}_{num}_{tend}.png",
                )

            # Create video
            if num == 100:
                video_path = results_dir / f"{xdmf_path.stem}_{num}_{tend}.mp4"
                gif_path = results_dir / f"{xdmf_path.stem}_{num}_{tend}.gif"
                create_video(
                    image_pattern=str(output_dir / "horizontal" / "sim_%05d.png"),
                    video_path=video_path,
                )
                create_gif_from_video(video_path=video_path, gif_path=gif_path)


if __name__ == "__main__":
    # -----------------------------------
    # Substrate scan
    # -----------------------------------
    # process files
    xdmf_dir = Path("/home/mkoenig/git/porous_media/data/spt_substrate_scan")
    xdmf_paths: List[Path] = xdmfs_from_febio(
        febio_dir=Path(
            "/home/mkoenig/git/porous_media/data/spt/simulation_fixedcelltype"
        ),
        xdmf_dir=xdmf_dir,
        overwrite=False,
    )
    info: XDMFInformation = XDMFInformation.from_path(xdmf_path=xdmf_paths[0])
    console.print(info)

    # create visualizations
    from porous_media.analyses.spt import data_layers_spt

    results_dir: Path = BASE_DIR / "results" / "spt_substrate_scan"
    visualize_scan(
        xdmf_paths=xdmf_paths,
        data_layers=data_layers_spt,
        results_dir=results_dir,
        create_panels=True,
    )

    # 5mmHg and 1mmHg => check effective_fluid_pressure!
    # perfusion = 1.0  /min
