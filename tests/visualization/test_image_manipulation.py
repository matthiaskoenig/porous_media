"""Check image composition with unequal panel dimensions."""

from pathlib import Path

import pytest
from PIL import Image

from porous_media.visualization.image_manipulation import merge_images


@pytest.fixture
def panels(tmp_path: Path) -> list[Path]:
    """Create differently sized panels with distinct colors."""
    paths = []
    for index, (size, color) in enumerate(
        [((4, 7), "red"), ((2, 3), "green"), ((3, 5), "blue")]
    ):
        path = tmp_path / f"panel{index}.png"
        with Image.new("RGB", size, color) as image:
            image.save(path)
        paths.append(path)
    return paths


@pytest.mark.parametrize(
    ("direction", "size", "green_origin", "blue_origin"),
    [
        ("vertical", (4, 15), (0, 7), (0, 10)),
        ("horizontal", (9, 7), (4, 0), (6, 0)),
        ("square", (8, 14), (4, 0), (0, 7)),
        ("custom", (12, 7), (4, 0), (8, 0)),
    ],
)
def test_merge_images_places_unequal_panels_without_overlap(
    tmp_path: Path,
    panels: list[Path],
    direction: str,
    size: tuple[int, int],
    green_origin: tuple[int, int],
    blue_origin: tuple[int, int],
) -> None:
    """Keep full panel heights and align grid panels with consistent spacing."""
    output = tmp_path / "merged.png"
    merge_images(iter(panels), output, direction=direction, ncols=3, nrows=1)
    with Image.open(output) as image:
        assert image.size == size
        assert image.getpixel((0, 0)) == (255, 0, 0)
        assert image.getpixel((3, 6)) == (255, 0, 0)
        assert image.getpixel(green_origin) == (0, 128, 0)
        assert image.getpixel(blue_origin) == (0, 0, 255)


@pytest.mark.parametrize(
    ("ncols", "nrows"), [(None, 1), (1, None), (0, 1), (1, -1), (1, 1)]
)
def test_merge_images_rejects_invalid_custom_grid(
    tmp_path: Path, panels: list[Path], ncols: int | None, nrows: int | None
) -> None:
    """Reject missing, nonpositive, or insufficient grid dimensions."""
    output = tmp_path / "invalid.png"
    with pytest.raises(ValueError):
        merge_images(panels, output, direction="custom", ncols=ncols, nrows=nrows)
    assert not output.exists()


def test_merge_images_rejects_empty_input(tmp_path: Path) -> None:
    """Explain an empty input without creating an image."""
    with pytest.raises(ValueError, match=r"[Aa]t least one image"):
        merge_images([], tmp_path / "empty.png")
