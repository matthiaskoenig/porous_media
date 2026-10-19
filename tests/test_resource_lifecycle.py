"""Exercise cleanup and bounded resource use through public workflows."""

import json
import subprocess
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pytest
from PIL import Image

from porous_media.log import get_logger
from porous_media.mesh.mesh_geojson import plot_geojson


def test_logger_lookup_does_not_duplicate_output(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """Repeatedly obtaining a named logger must emit one copy of each message."""
    name = f"review.{tmp_path.name}"
    logger = get_logger(name)
    try:
        get_logger(name).info("single message")
        assert capsys.readouterr().out.count("single message") == 1
    finally:
        for handler in list(logger.handlers):
            logger.removeHandler(handler)
            handler.close()


def test_geojson_failure_closes_created_figure(tmp_path: Path) -> None:
    """A malformed geometry must not leave an unused Matplotlib figure alive."""
    source = tmp_path / "invalid.geojson"
    source.write_text(
        json.dumps(
            {
                "type": "FeatureCollection",
                "features": [
                    {
                        "geometry": {
                            "type": "Polygon",
                            "coordinates": [[["invalid", 0], [0, 1], [1, 0]]],
                        }
                    }
                ],
            }
        ),
        encoding="utf-8",
    )
    before = set(plt.get_fignums())
    try:
        with pytest.raises(ValueError):
            plot_geojson(source)
        assert set(plt.get_fignums()) == before
    finally:
        for number in set(plt.get_fignums()) - before:
            plt.close(number)


def test_merge_many_images_under_file_limit(tmp_path: Path) -> None:
    """Merge more images than the process file limit without retaining open inputs."""
    if sys.platform == "win32":
        pytest.skip("POSIX file descriptor limits are unavailable")
    for number in range(100):
        with Image.new("RGB", (1, 1), "white") as image:
            image.save(tmp_path / f"frame_{number}.png")
    code = """
import resource
import sys
from pathlib import Path
from porous_media.visualization.image_manipulation import merge_images
root = Path(sys.argv[1])
_, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
resource.setrlimit(resource.RLIMIT_NOFILE, (32, hard))
merge_images(sorted(root.glob('frame_*.png')), root / 'merged.png')
"""
    result = subprocess.run(
        [sys.executable, "-c", code, str(tmp_path)],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    with Image.open(tmp_path / "merged.png") as merged:
        assert merged.size == (1, 100)
