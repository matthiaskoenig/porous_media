"""Image manipulation.

Combine and annotate images using PIL.
"""
from pathlib import Path
from typing import Iterable

from PIL import Image


def merge_images(paths: Iterable[Path], output_path: Path, direction: str = "vertical"):
    """Merge images either vertical or horizontal."""
    images = [Image.open(x) for x in paths]
    widths, heights = zip(*(i.size for i in images))

    if direction == "horizontal":
        total_width = sum(widths)
        max_height = max(heights)
        new_im = Image.new('RGB', (total_width, max_height))

        x_offset = 0
        for im in images:
          new_im.paste(im, (x_offset, 0))
          x_offset += im.size[0]

    elif direction == "vertical":
        total_height = sum(heights)
        max_width = max(widths)
        new_im = Image.new('RGB', (max_width, total_height))

        y_offset = 0
        for im in images:
            new_im.paste(im, (0, y_offset))
            y_offset += im.size[0]

    new_im.save(output_path)


if __name__ == "__main__":
    paths = ["test_pyvista_rr_(S).png", "test_pyvista_rr_(P).png", "test_pyvista_rr_necrosis.png"]
    merge_images(paths=paths, output_path="test_vertical.png", direction="vertical")
    merge_images(paths=paths, output_path="test_horizontal.png", direction="horizontal")
