"""Image manipulation.

Combine and annotate images using PIL.
"""
from pathlib import Path
from typing import Iterable, Optional, Tuple

import numpy as np
from PIL import Image

from porous_media.console import console


def merge_images(
    paths: Iterable[Path],
    output_path: Path,
    direction: str = "vertical",
    ncols: Optional[int] = None,
    nrows: Optional[int] = None,
) -> None:
    """Merge/combine images either vertical or horizontal or square.

    This creates larger images from individual panels. Panels can be combined
    vertical, horizonal, square or in custom layout.
    In case of custom layout height and width are required.
    """
    supported_directions = ["vertical", "horizontal", "square", "custom"]
    if direction not in supported_directions:
        raise ValueError(
            f"direction '{direction}' not in supported directions: {supported_directions}"
        )

    images = [Image.open(x) for x in paths]
    widths, heights = zip(*(i.size for i in images))

    if direction == "horizontal":
        total_width = sum(widths)
        max_height = max(heights)
        new_im = Image.new("RGB", (total_width, max_height))

        x_offset = 0
        for im in images:
            new_im.paste(im, (x_offset, 0))
            x_offset += im.size[0]

    elif direction == "vertical":
        total_height = sum(heights)
        max_width = max(widths)
        new_im = Image.new("RGB", (max_width, total_height))

        y_offset = 0
        for im in images:
            new_im.paste(im, (0, y_offset))
            y_offset += im.size[1]

    elif direction == "square":
        n = int(np.ceil(np.sqrt(len(images))))
        max_height = max(widths)
        max_width = max(widths)
        total_height = max_height * n
        total_width = max_width * n

        new_im = Image.new("RGB", (total_width, total_height))

        for k, im in enumerate(images):
            kx = int(k % n)
            ky = int(np.floor(k / n))
            x_offset = kx * im.size[0]
            y_offset = ky * im.size[1]
            new_im.paste(im, (x_offset, y_offset))

    elif direction == "custom":
        if not ncols:
            raise ValueError(f"width is required for direction '{direction}'")
        if not nrows:
            raise ValueError(f"height is required for direction '{direction}'")

        max_height = max(widths)
        max_width = max(widths)
        total_height = max_height * nrows
        total_width = max_width * ncols
        new_im = Image.new("RGB", (total_width, total_height))

        for k, im in enumerate(images):
            kx = int(k % ncols)
            ky = int(np.floor(k / ncols))
            x_offset = kx * im.size[0]
            y_offset = ky * im.size[1]
            new_im.paste(im, (x_offset, y_offset))

    new_im.save(output_path)


if __name__ == "__main__":
    # annotate image

    # merge images
    paths = [
        Path(p)
        for p in [
            "test_pyvista_rr_(S).png",
            "test_pyvista_rr_(P).png",
            "test_pyvista_rr_necrosis.png",
        ]
    ]
    merge_images(
        paths=paths, output_path=Path("test_vertical.png"), direction="vertical"
    )
    merge_images(
        paths=paths, output_path=Path("test_horizontal.png"), direction="horizontal"
    )
