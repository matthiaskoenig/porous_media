"""2D visualisation of results in lobulus geometry.

Creates static plots and videos of the TPM simulations.
"""

from pathlib import Path
from typing import Iterable

import numpy as np

from porous_media.data.xdmf_tools import (
    DataLimits,
    XDMFInfo,
    interpolate_xdmf,
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


def visualize_spt_2d(
    xdmf_paths: Iterable[Path],
    data_layers: list[DataLayer],
    results_dir: Path,
    selection: list[str],
    create_panels: bool = True,
) -> None:
    """Create static images and video."""

    # Calculate tend time from all simulations
    tends: np.ndarray = np.zeros(shape=(len(list(xdmf_paths)),))
    for k, xdmf_path in enumerate(xdmf_paths):
        xdmf_info = XDMFInfo.from_path(xdmf_path)
        tends[k] = xdmf_info.tend
    tend = tends.min()

    # ordered data layers selected
    data_layers_dict = {dl.sid: dl for dl in data_layers}
    data_layers_selected = [data_layers_dict[sid] for sid in selection]

    if create_panels:
        for num in [
            10,
            # 200
        ]:
            output_dir = results_dir / f"{num}_{tend}"
            output_dir.mkdir(exist_ok=True, parents=True)

            # interpolate and global limits
            all_limits: list[DataLimits] = []
            for xdmf_path in xdmf_paths:
                # interpolate & calculate limits
                interpolate_xdmf(
                    xdmf_in=xdmf_path,
                    xdmf_out=output_dir / f"{xdmf_path.stem}_interpolated.xdmf",
                    times_interpolate=np.linspace(0, tend, num=num),
                    overwrite=False,
                )

                limits = DataLimits.from_xdmf(xdmf_path=xdmf_path, overwrite=True)
                all_limits.append(limits)

            # merge limits from different simulations
            data_limits = DataLimits.merge_limits(all_limits)
            for data_layer in data_layers_selected:
                # only update empty limits
                data_layer.update_color_limits(data_limits=data_limits, only_empty=True)

            # create all panels for interpolation
            for xdmf_path in xdmf_paths:
                visualize_datalayers_timecourse(
                    xdmf_path=output_dir / f"{xdmf_path.stem}_interpolated.xdmf",
                    data_layers=data_layers_selected,
                    output_dir=output_dir / f"{xdmf_path.stem}",
                )

    # Create combined images for all simulations
    for num in [10, 40]:
        for xdmf_path in xdmf_paths:
            output_dir = results_dir / f"{num}_{tend}" / f"{xdmf_path.stem}"
            rows: list[Path] = create_combined_images(
                num_steps=num,
                direction="horizontal",
                output_dir=output_dir,
                selection=selection,
            )

            # Create combined figure for timecourse
            if num == 10:
                merge_images(
                    paths=[rows[k] for k in [1, 3, 5, 7, 9]],
                    direction="vertical",
                    output_path=results_dir / f"{xdmf_path.stem}_{num}_{tend}.png",
                )

            # Create video
            if num == 40:
                video_path = results_dir / f"{xdmf_path.stem}_{num}_{tend}.mp4"
                gif_path = results_dir / f"{xdmf_path.stem}_{num}_{tend}.gif"
                create_video(
                    image_pattern=str(output_dir / "horizontal" / "sim_%05d.png"),
                    video_path=video_path,
                )
                create_gif_from_video(video_path=video_path, gif_path=gif_path)


if __name__ == "__main__":
    from porous_media.analyses.spt import data_spt_dir

    analysis_folder = "convergence_sixth"
    # analysis_folder = "convergence_lobulus"

    xdmf_dir = data_spt_dir / analysis_folder / "xdmf"
    xdmf_paths = sorted([f for f in xdmf_dir.glob("*.xdmf")])

    results_dir: Path = data_spt_dir / analysis_folder / "2D"
    results_dir.mkdir(parents=True, exist_ok=True)

    # create visualizations
    from porous_media.analyses.spt import data_layers_spt, selection_spt

    visualize_spt_2d(
        xdmf_paths=xdmf_paths,
        data_layers=data_layers_spt,
        selection=selection_spt,
        results_dir=results_dir,
        create_panels=True,
    )
