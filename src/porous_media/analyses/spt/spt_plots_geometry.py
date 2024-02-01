"""2D visualisation of results in lobulus geometry.

Creates static plots and videos of the TPM simulations.
"""

from pathlib import Path
from typing import Dict, Iterable, List, Set

import numpy as np

from porous_media import BASE_DIR, DATA_DIR
from porous_media.console import console
from porous_media.data.xdmf_tools import (
    DataLimits,
    XDMFInfo,
    interpolate_xdmf,
    xdmfs_from_directory,
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
    data_layers: List[DataLayer],
    results_dir: Path,
    selection: List[str],
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
        for num in [10, 200]:  # [10, 100]:
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
    for num in [10, 200]:  # [10, 100]:
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
                    paths=[rows[k] for k in [1, 3, 5, 7, 9]],  # FIXME: hardcoded subset
                    direction="vertical",
                    output_path=results_dir / f"{xdmf_path.stem}_{num}_{tend}.png",
                )

            # Create video
            if num == 200:
                video_path = results_dir / f"{xdmf_path.stem}_{num}_{tend}.mp4"
                gif_path = results_dir / f"{xdmf_path.stem}_{num}_{tend}.gif"
                create_video(
                    image_pattern=str(output_dir / "horizontal" / "sim_%05d.png"),
                    video_path=video_path,
                )
                create_gif_from_video(video_path=video_path, gif_path=gif_path)


def necrosis_plots(xdmf_paths: List[Path], results_dir: Path) -> None:
    """Combine panels for the necrosis plots."""
    console.rule("Create necrosis plots", style="white", align="left")

    # Calculate tend time from all simulations
    tends: np.ndarray = np.zeros(shape=(len(list(xdmf_paths)),))
    for k, xdmf_path in enumerate(xdmf_paths):
        xdmf_info = XDMFInfo.from_path(xdmf_path)
        tends[k] = xdmf_info.tend
    tend = tends.min()

    # get all necrosis pictures
    num_time = 200
    num_substrate = 8
    num_patterns = 6

    necrosis_paths: List[Path] = []
    for k, xdmf_path in enumerate(xdmf_paths):
        p = (
            results_dir
            / f"{num_time}_{tend}"
            / f"{xdmf_path.stem}"
            / "panels"
            / "rr_necrosis"
            / f"sim_{num_time-1:05d}.png"
        )
        necrosis_paths.append(p)
        # add the zonation pattern at the end
        if (k > 0) & ((k + 1) % num_substrate == 0):
            p_pattern = (
                results_dir
                / f"{num_time}_{tend}"
                / f"{xdmf_path.stem}"
                / "panels"
                / "rr_protein"
                / f"sim_{num_time - 1:05d}.png"
            )
            necrosis_paths.append(p_pattern)

    image_path: Path = results_dir / "zonation_pattern_necrosis.png"
    merge_images(
        paths=necrosis_paths,
        direction="custom",
        ncols=num_substrate + 1,
        nrows=num_patterns,
        output_path=image_path,
    )
    console.print(f"file://{image_path}")


if __name__ == "__main__":
    # TODO: create videos for plots (40 x)
    # TODO: create timecourse plots (40 x) -> select one plot

    # TODO: combine static plots at 10 hr from panels

    # date = "2023-12-13"
    date = "2023-12-19"
    xdmf_dir = Path(f"/home/mkoenig/git/porous_media/data/spt/{date}/xdmf")
    xdmf_paths = sorted([f for f in xdmf_dir.glob("*.xdmf")])

    # create visualizations
    from porous_media.analyses.spt import data_layers_spt, selection_spt

    results_dir: Path = BASE_DIR / "results" / "spt" / date / "2D"
    results_dir.mkdir(parents=True, exist_ok=True)
    visualize_spt_2d(
        xdmf_paths=xdmf_paths,
        data_layers=data_layers_spt,
        selection=selection_spt,
        results_dir=results_dir,
        create_panels=True,
    )

    necrosis_plots(xdmf_paths=xdmf_paths, results_dir=results_dir)
