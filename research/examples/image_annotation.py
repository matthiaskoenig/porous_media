"""Archived example; configure research paths before running."""

from pathlib import Path

from porous_media import RESOURCES_DIR
from porous_media.visualization.image_annotation import annotate_image_text
from research import RESULTS_DIR

if __name__ == "__main__":
    image_in: Path = RESOURCES_DIR / "images" / "mesh_zonation.png"
    image_out: Path = RESULTS_DIR / "mesh_zonation_annotated.png"

    annotate_image_text(
        image_in_path=image_in, image_out_path=image_out, text="time = 123.2 [min]"
    )
