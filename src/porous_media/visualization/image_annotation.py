"""Annotate images with text."""
from pathlib import Path
from typing import Tuple

from PIL import Image, ImageDraw, ImageFont

from porous_media import RESOURCES_DIR, RESULTS_DIR
from porous_media.console import console


def annotate_image_text(
    image_in_path: Path, image_out_path: Path, text: str, xy: Tuple[int, int] = (10, 10)
) -> None:
    """Annotate images with text and store the resulting image."""
    console.print(f"Annotate: {image_in_path}")

    image = Image.open(image_in_path)
    draw = ImageDraw.Draw(image)
    font = ImageFont.truetype("arial.ttf", 40, encoding="unic")
    # console.print(font.get_variation_names())
    # font.set_variation_by_name('Bold')
    draw.text(xy, text, fill="#000000", font=font)
    image.save(image_out_path, "PNG")
    console.print(f"Annotated image: file://{image_out_path}")


if __name__ == "__main__":
    image_in: Path = RESOURCES_DIR / "images" / "mesh_zonation.png"
    image_out: Path = RESULTS_DIR / "mesh_zonation_annotated.png"

    annotate_image_text(
        image_in_path=image_in, image_out_path=image_out, text="time = 123.2 [min]"
    )
