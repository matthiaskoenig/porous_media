"""Archived example; configure research paths before running."""

from pathlib import Path

from porous_media.visualization.image_manipulation import merge_images

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
