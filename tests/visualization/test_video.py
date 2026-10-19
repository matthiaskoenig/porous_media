"""Exercise video helpers through the ffmpeg process boundary."""

import shutil
import subprocess
from pathlib import Path

import pytest
from PIL import Image

from porous_media.visualization.video import create_gif_from_video, create_video


@pytest.fixture
def video_file(tmp_path: Path) -> Path:
    """Encode a short video when the optional ffmpeg program is available."""
    if shutil.which("ffmpeg") is None:
        pytest.skip("ffmpeg is required for video integration tests")
    for index, color in enumerate(("red", "blue")):
        with Image.new("RGB", (16, 16), color) as image:
            image.save(tmp_path / f"frame{index}.png")
    video = tmp_path / "video output.mp4"
    create_video(str(tmp_path / "frame%d.png"), video, frame_rate=1, overwrite=False)
    assert video.stat().st_size > 0
    return video


def test_create_video_refuses_to_overwrite_existing_output(video_file: Path) -> None:
    """Reject overwriting and preserve the existing video."""
    original = video_file.read_bytes()
    with pytest.raises(FileExistsError):
        create_video(
            str(video_file.parent / "frame%d.png"), video_file, overwrite=False
        )
    assert video_file.read_bytes() == original


def test_create_gif_refuses_to_overwrite_existing_output(video_file: Path) -> None:
    """Create a GIF and propagate a refused overwrite without damaging the file."""
    gif = video_file.with_suffix(".gif")
    create_gif_from_video(video_file, gif, overwrite=False)
    original = gif.read_bytes()
    assert original.startswith(b"GIF")
    with pytest.raises(FileExistsError):
        create_gif_from_video(video_file, gif, overwrite=False)
    assert gif.read_bytes() == original


def test_create_video_reports_encoding_failure(video_file: Path) -> None:
    """Propagate ffmpeg errors when the image pattern matches no frames."""
    output = video_file.parent / "missing.mp4"
    with pytest.raises(subprocess.CalledProcessError):
        create_video(str(video_file.parent / "missing%d.png"), output)
    assert not output.exists()


def test_create_gif_rejects_missing_video(tmp_path: Path) -> None:
    """Report a missing input before invoking ffmpeg."""
    with pytest.raises(OSError, match="Video path does not exist"):
        create_gif_from_video(tmp_path / "missing.mp4", tmp_path / "output.gif")
