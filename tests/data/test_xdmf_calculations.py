"""Check scalar datasets with mixed cell blocks and reusable disk caches."""

from pathlib import Path

import meshio
import numpy as np
import pytest
import xarray as xr

from porous_media.data.xdmf_calculations import (
    mesh_datasets_from_xdmf,
    mesh_datasets_from_xdmf_dir,
)


def write_series(path: Path, *, changing_schema: bool = False) -> None:
    """Write unequal cell blocks and values that promote from integers to floats."""
    points = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0]]
    )
    cells = [("triangle", [[0, 1, 2], [1, 3, 2]]), ("quad", [[0, 1, 3, 2]])]
    with meshio.xdmf.TimeSeriesWriter(path, data_format="XML") as writer:
        writer.write_points_cells(points, cells)
        for step in range(2):
            values = (
                np.array([1, 2, 3], dtype=np.int32)
                if not step
                else np.array([1.25, 2.25, 3.25])
            )
            field = "changed" if step and changing_schema else "pressure"
            writer.write_data(
                float(step),
                point_data={
                    "temperature": np.arange(4) + step,
                    "vector": np.ones((4, 3)),
                },
                cell_data={
                    field: [values[:2], values[2:]],
                    "vector": [np.ones((2, 3)), np.ones((1, 3))],
                },
            )


def test_scalar_datasets_keep_all_blocks_and_fractional_values(tmp_path: Path) -> None:
    """Preserve mixed-block ordering, scalar dimensions, and numeric promotion."""
    source = tmp_path / "source.xdmf"
    write_series(source)
    cells, points = mesh_datasets_from_xdmf(source)
    assert cells["pressure"].dims == ("time", "cell")
    assert points["temperature"].dims == ("time", "point")
    np.testing.assert_allclose(cells["pressure"], [[1, 2, 3], [1.25, 2.25, 3.25]])
    np.testing.assert_array_equal(points["temperature"], [[0, 1, 2, 3], [1, 2, 3, 4]])
    assert "vector" not in cells and "vector" not in points
    cached_cells, cached_points = mesh_datasets_from_xdmf(source)
    xr.testing.assert_equal(cells, cached_cells)
    xr.testing.assert_equal(points, cached_points)
    source.with_name("source_cells.nc").unlink()
    source.with_name("source_points.nc").unlink()
    np.testing.assert_array_equal(cached_cells["pressure"][0], [1, 2, 3])


def test_scalar_datasets_keep_single_cell_and_empty_point_fields(
    tmp_path: Path,
) -> None:
    """Do not squeeze away a scalar cell or invent fields for empty point data."""
    source = tmp_path / "single.xdmf"
    with meshio.xdmf.TimeSeriesWriter(source, data_format="XML") as writer:
        writer.write_points_cells(np.eye(3), [("triangle", [[0, 1, 2]])])
        writer.write_data(0.0, cell_data={"value": [np.array([7.0])]})
    cells, points = mesh_datasets_from_xdmf(source)
    assert cells["value"].shape == (1, 1)
    assert points.sizes == {"time": 1, "point": 3}
    assert not points.data_vars


def test_scalar_datasets_reject_changing_fields_without_caches(tmp_path: Path) -> None:
    """Report inconsistent inputs before writing a misleading cache pair."""
    source = tmp_path / "source.xdmf"
    write_series(source, changing_schema=True)
    with pytest.raises(ValueError, match="consistent"):
        mesh_datasets_from_xdmf(source)
    assert not list(tmp_path.glob("*.nc"))


@pytest.mark.parametrize("workers", [1, 2])
def test_directory_conversion_builds_caches(tmp_path: Path, workers: int) -> None:
    """Produce equivalent files with serial or bounded process execution."""
    for name in ("first", "second"):
        write_series(tmp_path / f"{name}.xdmf")
    mesh_datasets_from_xdmf_dir(tmp_path, max_workers=workers)
    for name in ("first", "second"):
        assert (tmp_path / f"{name}_cells.nc").is_file()
        assert (tmp_path / f"{name}_points.nc").is_file()


def test_directory_conversion_validates_workers(tmp_path: Path) -> None:
    """Allow an empty directory while rejecting invalid worker limits."""
    mesh_datasets_from_xdmf_dir(tmp_path)
    with pytest.raises(ValueError, match="max_workers"):
        mesh_datasets_from_xdmf_dir(tmp_path, max_workers=0)
