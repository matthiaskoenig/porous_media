"""Tools related to videos.

Requires `ffmpeg` for video generation.
"""

import os
from pathlib import Path

from porous_media.console import console


def create_video(
    image_pattern: str, video_path: Path, frame_rate: int = 30, codec: str = "mpeg4"
) -> None:
    """Create video by combining images."""
    command = f"ffmpeg -f image2 -r {frame_rate} -i {image_pattern} -vcodec {codec} -y {video_path}"
    console.print(f"Create video: {video_path}")
    console.print(command)
    os.system(command)


def create_gif_from_video(video_path: Path, gif_path: Path) -> None:
    """Gif from video."""
    if not video_path.exists():
        raise IOError(f"Video path does not exist: {video_path}")

    command = f"ffmpeg -i {video_path} {gif_path}"
    console.print(f"Create gif: {gif_path}")
    console.print(command)
    os.system(command)


if __name__ == "__main__":
    from porous_media import RESULTS_DIR

    video_path = (
        RESULTS_DIR / "spt_zonation_patterns" / "simulation_pattern0_100_28800.mp4"
    )
    gif_path = (
        RESULTS_DIR / "spt_zonation_patterns" / "simulation_pattern0_100_28800.gif"
    )
    create_gif_from_video(
        video_path=video_path,
        gif_path=gif_path,
    )
