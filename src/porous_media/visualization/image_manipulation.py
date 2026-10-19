"""Image manipulation.

Combine and annotate images using PIL.
"""

from collections.abc import Iterable
from math import ceil, sqrt
from pathlib import Path

from PIL import Image


def merge_images(
    paths: Iterable[Path],
    output_path: Path,
    direction: str = "vertical",
    ncols: int | None = None,
    nrows: int | None = None,
) -> None:
    """Combine images into a vertical strip, horizontal strip, or grid.

    Args:
        paths: Input images in panel order.
        output_path: Destination image path.
        direction: Layout: "vertical", "horizontal", "square", or "custom".
        ncols: Positive column count required for a custom grid.
        nrows: Positive row count required for a custom grid.

    Raises:
        ValueError: If no images are supplied or the layout is invalid.
    """
    supported_directions = ["vertical", "horizontal", "square", "custom"]
    if direction not in supported_directions:
        raise ValueError(
            f"direction '{direction}' not in supported directions: {supported_directions}"
        )

    paths = tuple(paths)
    if not paths:
        raise ValueError("paths must contain at least one image")

    sizes = []
    for path in paths:
        with Image.open(path) as image:
            sizes.append(image.size)
    widths, heights = zip(*sizes, strict=True)
    max_width, max_height = max(widths), max(heights)
    offsets: list[tuple[int, int]] = []
    if direction == "horizontal":
        size = (sum(widths), max_height)
        offset = 0
        for width in widths:
            offsets.append((offset, 0))
            offset += width
    elif direction == "vertical":
        size = (max_width, sum(heights))
        offset = 0
        for height in heights:
            offsets.append((0, offset))
            offset += height
    else:
        if direction == "square":
            ncols = nrows = ceil(sqrt(len(paths)))
        if ncols is None or ncols <= 0:
            raise ValueError("ncols must be positive for a custom grid")
        if nrows is None or nrows <= 0:
            raise ValueError("nrows must be positive for a custom grid")
        if ncols * nrows < len(paths):
            raise ValueError("custom grid must have room for every image")

        size = (max_width * ncols, max_height * nrows)
        offsets = [
            ((k % ncols) * max_width, (k // ncols) * max_height)
            for k in range(len(paths))
        ]

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with Image.new("RGB", size) as merged:
        for path, position in zip(paths, offsets, strict=True):
            with Image.open(path) as image:
                merged.paste(image, position)
        merged.save(output_path)
