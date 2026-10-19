"""Annotate images with text."""

from pathlib import Path

from PIL import Image, ImageDraw, ImageFont

from porous_media.console import console


def annotate_image_text(
    image_in_path: Path,
    image_out_path: Path,
    text: str,
    xy: tuple[int, int] = (10, 10),
    *,
    font_path: Path | None = None,
    font_size: int = 40,
) -> None:
    """Annotate an image with a portable default font or an explicit TrueType font."""
    if font_size <= 0:
        raise ValueError("font_size must be positive")
    console.print(f"Annotate: {image_in_path}")

    font = (
        ImageFont.load_default(size=font_size)
        if font_path is None
        else ImageFont.truetype(str(font_path), font_size)
    )
    with Image.open(image_in_path) as image:
        draw = ImageDraw.Draw(image)
        draw.text(xy, text, fill="#000000", font=font)
        image_out_path.parent.mkdir(parents=True, exist_ok=True)
        image.save(image_out_path, "PNG")
    console.print(f"Annotated image: file://{image_out_path}")
