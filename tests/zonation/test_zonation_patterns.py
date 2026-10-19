"""Check spatial correctness and normalization independently of graphics."""

from pathlib import Path

import meshio
import numpy as np
import pytest

from porous_media.mesh.mesh_zonation import ZonatedMesh, ZonationPatterns


def boundary_mesh() -> meshio.Mesh:
    """Create mixed blocks whose boundary centers are ten units apart."""
    points = np.array(
        [
            [-0.1, 0.0, 0.0],
            [0.1, 0.0, 0.0],
            [1.9, 0.0, 0.0],
            [2.1, 0.0, 0.0],
            [9.9, 0.0, 0.0],
            [10.1, 0.0, 0.0],
        ]
    )
    return meshio.Mesh(
        points,
        [("line", [[0, 1], [2, 3]]), ("line", [[4, 5]])],
        cell_data={
            "cell_type": [np.array([1, 0]), np.array([2])],
            "element_volume_point_TPM": [np.array([1.0, 2.0]), np.array([3.0])],
        },
    )


@pytest.mark.parametrize("scale", [1.0, 1000.0])
def test_zonation_uses_all_blocks_and_true_distances(scale: float) -> None:
    """Preserve block boundaries and avoid clipping distances to one mesh unit."""
    mesh = boundary_mesh()
    mesh.points *= scale
    result = ZonatedMesh.create_zonated_mesh(mesh)
    np.testing.assert_allclose(np.asarray(result.cell_data["position"][0]), [0.0, 0.2])
    np.testing.assert_allclose(np.asarray(result.cell_data["position"][1]), [1.0])
    assert "position" not in mesh.cell_data


def test_random_pattern_is_reproducible_and_normalized(tmp_path: Path) -> None:
    """Normalize the actual random sample rather than a second independent sample."""
    source = tmp_path / "mesh.vtk"
    meshio.write(source, boundary_mesh())
    results = [
        ZonatedMesh.create_zonated_mesh_from_vtk(
            source, patterns=[ZonationPatterns.random], rng=np.random.default_rng(42)
        )
        for _ in range(2)
    ]
    values = np.concatenate(results[0].cell_data["pattern__random"])
    volume = np.concatenate(results[0].cell_data["volume"])
    assert np.average(values, weights=volume) == pytest.approx(0.5)
    np.testing.assert_array_equal(
        values, np.concatenate(results[1].cell_data["pattern__random"])
    )


@pytest.mark.parametrize("bad_volume", [0.0, -1.0, np.nan, np.inf])
def test_zonation_rejects_invalid_volumes_before_mutation(bad_volume: float) -> None:
    """Leave an input mesh unchanged when normalization cannot be defined."""
    mesh = boundary_mesh()
    mesh.cell_data["element_volume_point_TPM"] = [
        np.full(2, bad_volume),
        np.full(1, bad_volume),
    ]
    with pytest.raises(ValueError, match="volumes"):
        ZonatedMesh.create_zonated_mesh(mesh, copy_mesh=False)
    assert "element_volume_point_TPM" in mesh.cell_data
    assert "position" not in mesh.cell_data


def test_zonation_requires_both_boundaries() -> None:
    """Reject incomplete labels rather than inventing a boundary distance."""
    mesh = boundary_mesh()
    mesh.cell_data["cell_type"][1] = np.zeros(1)
    with pytest.raises(ValueError, match="Both periportal and perivenous"):
        ZonatedMesh.create_zonated_mesh(mesh)


def test_preserved_fields_accept_original_volume_name() -> None:
    """Use the documented original volume field even when preserving other fields."""
    mesh = boundary_mesh()
    result = ZonatedMesh.create_zonated_mesh(mesh, remove_cell_data=False)
    assert "element_volume_point_TPM" in result.cell_data
    np.testing.assert_array_equal(result.cell_data["volume"][0], [1.0, 2.0])


def test_exponential_decrease_respects_requested_endpoints() -> None:
    """Handle ranges whose minimum and maximum do not sum to one."""
    values = ZonationPatterns.exp_decrease(
        np.array([0.0, 1.0]), value_min=0.2, value_max=0.6
    )
    np.testing.assert_allclose(values, [0.6, 0.2])
