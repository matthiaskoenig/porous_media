"""Check point deduplication and preservation of mesh data."""

from pathlib import Path

import meshio
import numpy as np
import pytest

from porous_media.analyses.geometry_conversion import (
    reconstruct_lobulus_from_hexagon,
    remove_duplicate_points,
    unique_with_tolerance,
)


def test_unique_with_tolerance_groups_first_duplicate_once() -> None:
    """Map identical first points to one group instead of treating zero as unassigned."""
    points = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    unique, inverse = unique_with_tolerance(points, tol=1e-6)
    np.testing.assert_allclose(unique, points[[0, 2]])
    np.testing.assert_array_equal(inverse, [0, 0, 1])


def test_remove_duplicate_points_preserves_scalar_and_vector_data() -> None:
    """Select original data by each group's first point and remap connectivity."""
    points = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    mesh = meshio.Mesh(
        points,
        [("line", np.array([[0, 2], [1, 2]]))],
        point_data={
            "scalar": np.array([10, 20, 30]),
            "vector": np.array([[1, 2], [3, 4], [5, 6]]),
        },
        cell_data={"region": [np.array([7, 8])]},
    )
    result = remove_duplicate_points(mesh, points[[0, 2]], np.array([0, 0, 1]))
    np.testing.assert_array_equal(result.cells[0].data, [[0, 1], [0, 1]])
    np.testing.assert_array_equal(result.point_data["scalar"], [10, 30])
    np.testing.assert_array_equal(result.point_data["vector"], [[1, 2], [5, 6]])
    np.testing.assert_array_equal(result.cell_data["region"][0], [7, 8])
    np.testing.assert_array_equal(mesh.point_data["scalar"], [10, 20, 30])


@pytest.mark.parametrize("output_directory", [".", "nested"])
def test_reconstruct_lobulus_preserves_blocks_data_and_hdf_location(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, output_directory: str
) -> None:
    """Roundtrip repeated cell types and scalar/vector data without shared HDF files."""
    monkeypatch.chdir(tmp_path)
    source = tmp_path / "source.xdmf"
    points = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0]]
    )
    cells = [
        ("triangle", np.array([[0, 1, 2], [1, 3, 2]])),
        ("quad", np.array([[0, 1, 3, 2]])),
        ("triangle", np.array([[0, 2, 3]])),
    ]
    with meshio.xdmf.TimeSeriesWriter(source, data_format="XML") as writer:
        writer.write_points_cells(points, cells)
        writer.write_data(
            1.5,
            point_data={
                "scalar": np.array([1.0, 2.0, 3.0, 4.0]),
                "vector": np.arange(12).reshape(4, 3),
            },
            cell_data={
                "region": [np.array([10, 11]), np.array([20]), np.array([30])],
                "direction": [
                    np.ones((2, 3)),
                    2 * np.ones((1, 3)),
                    3 * np.ones((1, 3)),
                ],
            },
        )

    output = Path(output_directory) / "reconstructed.xdmf"
    unrelated_hdf = tmp_path / "reconstructed.h5"
    if output_directory != ".":
        unrelated_hdf.write_bytes(b"unrelated simulation")
    reconstruct_lobulus_from_hexagon(source, output)

    assert output.with_suffix(".h5").is_file()
    if output_directory != ".":
        assert unrelated_hdf.read_bytes() == b"unrelated simulation"
    with meshio.xdmf.TimeSeriesReader(output) as reader:
        result_points, result_cells = reader.read_points_cells()
        time, point_data, cell_data = reader.read_data(0)

    assert time == pytest.approx(1.5)
    assert result_points.shape == (24, 3)
    assert [block.type for block in result_cells] == ["triangle", "quad", "triangle"]
    assert [len(block) for block in result_cells] == [12, 6, 6]
    for original, result in zip(cells, result_cells, strict=True):
        original_data = original[1]
        np.testing.assert_array_equal(result.data[: len(original_data)], original_data)
        np.testing.assert_array_equal(
            result.data[-len(original_data) :], original_data + 20
        )
    np.testing.assert_array_equal(point_data["scalar"], [1.0, 2.0, 3.0, 4.0] * 6)
    np.testing.assert_array_equal(
        point_data["vector"], np.tile(np.arange(12).reshape(4, 3), (6, 1))
    )
    np.testing.assert_array_equal(cell_data["region"][0], [10, 11] * 6)
    np.testing.assert_array_equal(cell_data["region"][1], [20] * 6)
    np.testing.assert_array_equal(cell_data["region"][2], [30] * 6)
    for index, block in enumerate(cell_data["direction"]):
        np.testing.assert_array_equal(
            block, np.full((len(result_cells[index]), 3), index + 1)
        )
