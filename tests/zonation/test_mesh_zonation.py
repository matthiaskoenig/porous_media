"""Test that all resources for the app are available."""
from pm.mesh.mesh_zonation import example_mesh_zonation


def test_mesh_zonation_example(tmp_path) -> None:
    """Test example."""
    example_mesh_zonation(results_dir=tmp_path)
