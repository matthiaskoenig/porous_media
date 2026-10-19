"""Archived example; configure research paths before running."""

from porous_media.mesh.mesh_zonation import example_mesh_zonation

if __name__ == "__main__":
    from research import RESULTS_DIR

    example_mesh_zonation(results_dir=RESULTS_DIR)
