"""Test that all resources for the app are available."""

from pathlib import Path

import meshio
import numpy as np
import pytest

from porous_media.mesh import mesh_zonation
from porous_media.mesh.mesh_zonation import example_mesh_zonation


def test_mesh_zonation_example(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """The nonvisual example exports zonation data without opening a plotter."""

    def unexpected_visualization(*args: object, **kwargs: object) -> None:
        pytest.fail("visualize=False must not start visualization")

    monkeypatch.setattr(
        mesh_zonation, "visualize_zonation_patterns", unexpected_visualization
    )
    example_mesh_zonation(results_dir=tmp_path, visualize=False)

    output_dir = tmp_path / "mesh_zonation"
    assert (output_dir / "mesh_zonation_lobulus.xdmf").is_file()
    mesh = meshio.read(output_dir / "mesh_zonation_lobulus.vtk")
    position = mesh.cell_data["position"][0]
    assert np.isfinite(position).all()
    assert ((position >= 0.0) & (position <= 1.0)).all()
