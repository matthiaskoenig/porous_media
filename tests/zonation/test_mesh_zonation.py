"""Test that all resources for the app are available."""
from pathlib import Path

from porous_media.mesh.mesh_zonation import example_mesh_zonation


def test_mesh_zonation_example(tmp_path: Path) -> None:
    """Test example."""
    example_mesh_zonation(results_dir=tmp_path, visualize=False)
