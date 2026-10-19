"""Tools related to videos.

Requires `ffmpeg` for video generation.
"""

import subprocess
from pathlib import Path

from porous_media.console import console


def create_video(
    image_pattern: str,
    video_path: Path,
    frame_rate: int = 30,
    codec: str = "mpeg4",
    overwrite: bool = True,
) -> None:
    """Create a video from image frames using ffmpeg.

    Args:
        image_pattern: ffmpeg input pattern for the image sequence.
        video_path: Destination video path.
        frame_rate: Frames per second.
        codec: ffmpeg video codec name.
        overwrite: Whether ffmpeg may overwrite an existing video.

    Raises:
        FileExistsError: If the destination exists and overwrite is false.
        subprocess.CalledProcessError: If ffmpeg fails.
    """
    if not overwrite and video_path.exists():
        raise FileExistsError(f"Video path already exists: {video_path}")

    overwrite_flag = "-y" if overwrite else "-n"
    command = [
        "ffmpeg",
        overwrite_flag,
        "-f",
        "image2",
        "-r",
        str(frame_rate),
        "-i",
        image_pattern,
        "-vcodec",
        codec,
        video_path,
    ]
    console.print(f"Create video: {video_path}")
    console.print(command)
    subprocess.run(command, check=True)


def create_gif_from_video(
    video_path: Path, gif_path: Path, overwrite: bool = True
) -> None:
    """Convert a video to a GIF using ffmpeg.

    Args:
        video_path: Existing source video path.
        gif_path: Destination GIF path.
        overwrite: Whether ffmpeg may overwrite an existing GIF.

    Raises:
        OSError: If the source video does not exist.
        FileExistsError: If the destination exists and overwrite is false.
        subprocess.CalledProcessError: If ffmpeg fails.
    """
    if not video_path.exists():
        raise OSError(f"Video path does not exist: {video_path}")
    if not overwrite and gif_path.exists():
        raise FileExistsError(f"GIF path already exists: {gif_path}")

    overwrite_flag = "-y" if overwrite else "-n"
    command = ["ffmpeg", overwrite_flag, "-i", video_path, gif_path]
    console.print(f"Create gif: {gif_path}")
    console.print(command)
    subprocess.run(command, check=True)


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
