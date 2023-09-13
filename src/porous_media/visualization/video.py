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
