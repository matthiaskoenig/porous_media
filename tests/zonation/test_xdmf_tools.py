"""Test xdmf functionality."""

from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import meshio
import numpy as np
import pytest

from porous_media import RESOURCES_DIR
from porous_media.data.xdmf_tools import (
    XDMFInfo,
    XDMFTimeSeriesWriter,
    interpolate_xdmf,
    vtks_to_xdmf,
)


def _write_time_series(path: Path, times: tuple[float, ...]) -> None:
    """Write a small portable time series with differently sized cell blocks."""
    points = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0]]
    )
    cells = [
        ("triangle", np.array([[0, 1, 2], [1, 3, 2]])),
        ("quad", np.array([[0, 1, 3, 2]])),
    ]
    with meshio.xdmf.TimeSeriesWriter(path, data_format="XML") as writer:
        writer.write_points_cells(points, cells)
        for time in times:
            writer.write_data(
                time,
                point_data={"pressure": np.full(4, 10 * time)},
                cell_data={
                    "concentration": [
                        np.full(2, 10 + 10 * time),
                        np.array([20 + 20 * time]),
                    ]
                },
            )


@pytest.fixture
def time_series(tmp_path: Path) -> Path:
    """Provide a time series that needs no HDF sidecar."""
    source = tmp_path / "source.xdmf"
    _write_time_series(source, (0.0, 10.0))
    return source


def test_xdmf_writer_produces_relocatable_hdf(tmp_path: Path) -> None:
    """Keep meshio's context manager API and relative references to the HDF sidecar."""
    original = tmp_path / "original"
    original.mkdir()
    output = original / "simulation.xdmf"
    with XDMFTimeSeriesWriter(output) as writer:
        writer.write_points_cells(
            np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
            [("triangle", np.array([[0, 1, 2]]))],
        )
        writer.write_data(0.0, point_data={"pressure": np.array([1.0, 2.0, 3.0])})

    relocated = tmp_path / "relocated"
    original.rename(relocated)
    with meshio.xdmf.TimeSeriesReader(relocated / output.name) as reader:
        reader.read_points_cells()
        time, point_data, _ = reader.read_data(0)
    assert time == pytest.approx(0.0)
    np.testing.assert_array_equal(point_data["pressure"], [1.0, 2.0, 3.0])


def test_vtk_single_to_xdmf(tmp_path: Path) -> None:
    """Test the vtk time course parsing."""
    vtk_dir = RESOURCES_DIR / "vtk" / "vtk_single"
    xdmf_path = tmp_path / "vtk_single.xdmf"
    vtks_to_xdmf(vtk_dir, xdmf_path=xdmf_path, overwrite=True)
    assert xdmf_path.exists()

    xdmf_info: XDMFInfo = XDMFInfo.from_path(xdmf_path)
    assert xdmf_info.num_steps == 1
    assert xdmf_info.tstart == pytest.approx(220.0)
    assert xdmf_info.tend == pytest.approx(220.0)

    for key in [
        "displacement",
        "effective_fluid_pressure_TPM",
    ]:
        assert key in xdmf_info.point_data

    for key in [
        "stress",
        "fluid_flux_TPM",
        "solid_stress_TPM",
        "S_ext",
        "P_ext",
    ]:
        assert key in xdmf_info.cell_data


def test_vtk_timecourse_to_xdmf(tmp_path: Path) -> None:
    """Test the vtk time course parsing."""
    vtk_dir = RESOURCES_DIR / "vtk" / "vtk_timecourse"
    xdmf_path = tmp_path / "vtk_timecourse.xdmf"
    vtks_to_xdmf(vtk_dir, xdmf_path=xdmf_path, overwrite=True)
    assert xdmf_path.exists()

    xdmf_info: XDMFInfo = XDMFInfo.from_path(xdmf_path)
    assert xdmf_info.num_steps == 3
    assert xdmf_info.tstart == pytest.approx(0.0)
    assert xdmf_info.tend == pytest.approx(220.0)


@pytest.mark.parametrize("output_directory", [".", "nested"])
def test_vtk_to_xdmf_keeps_hdf_next_to_output(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, output_directory: str
) -> None:
    """Conversion works in the current directory and preserves unrelated HDF files."""
    monkeypatch.chdir(tmp_path)
    output_dir = Path(output_directory)
    output_path = output_dir / "simulation.xdmf"
    unrelated_hdf = tmp_path / "simulation.h5"
    if output_directory != ".":
        unrelated_hdf.write_bytes(b"unrelated simulation")

    vtks_to_xdmf(RESOURCES_DIR / "vtk" / "vtk_single", output_path)

    assert output_path.with_suffix(".h5").is_file()
    assert XDMFInfo.from_path(output_path).num_steps == 1
    if output_directory != ".":
        assert unrelated_hdf.read_bytes() == b"unrelated simulation"


def test_interpolate_xdmf_interpolates_points_and_cell_blocks(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, time_series: Path
) -> None:
    """Interpolation preserves endpoints and linearly interpolates every cell block."""
    monkeypatch.chdir(tmp_path)
    output = tmp_path / "nested" / "interpolated.xdmf"
    interpolate_xdmf(time_series, output, np.array([0.0, 2.0, 10.0]))

    assert output.with_suffix(".h5").is_file()
    assert not (tmp_path / "interpolated.h5").exists()
    assert XDMFInfo.from_path(output).num_cells == 3
    with meshio.xdmf.TimeSeriesReader(output) as reader:
        reader.read_points_cells()
        for step, expected_time in enumerate((0.0, 2.0, 10.0)):
            time, point_data, cell_data = reader.read_data(step)
            assert time == pytest.approx(expected_time)
            np.testing.assert_allclose(point_data["pressure"], 10 * expected_time)
            np.testing.assert_allclose(
                cell_data["concentration"][0], 10 + 10 * expected_time
            )
            np.testing.assert_allclose(
                cell_data["concentration"][1], 20 + 20 * expected_time
            )


@pytest.mark.parametrize(
    ("times", "message"),
    [
        ([], "nonempty one-dimensional"),
        ([[1.0]], "nonempty one-dimensional"),
        ([np.nan], "must be finite"),
        ([np.inf], "must be finite"),
        ([0.0, -1.0, 10.0], "Lower interpolation range"),
        ([0.0, 11.0, 10.0], "Upper interpolation range"),
    ],
)
def test_interpolate_xdmf_rejects_invalid_times_without_output(
    tmp_path: Path, time_series: Path, times: list, message: str
) -> None:
    """Reject invalid times before creating any output files."""
    output = tmp_path / "nested" / "invalid.xdmf"
    with pytest.raises(ValueError, match=message):
        interpolate_xdmf(time_series, output, np.array(times))
    assert not output.parent.exists()


@pytest.mark.parametrize("times", [(0.0, 0.0), (10.0, 0.0), (0.0, np.nan)])
def test_interpolate_xdmf_rejects_invalid_source_times(
    tmp_path: Path, times: tuple[float, ...]
) -> None:
    """Require finite strictly increasing source times for unambiguous interpolation."""
    source = tmp_path / "source.xdmf"
    _write_time_series(source, times)
    output = tmp_path / "invalid.xdmf"
    with pytest.raises(ValueError, match="Input timepoints"):
        interpolate_xdmf(source, output, np.array([0.0]))
    assert not output.exists()
    assert not output.with_suffix(".h5").exists()


@pytest.mark.parametrize("times", [(0.0, 1e-10), (3.0,)])
def test_interpolate_xdmf_handles_small_intervals_and_single_timepoint(
    tmp_path: Path, times: tuple[float, ...]
) -> None:
    """Preserve interpolation precision for short intervals and single-step inputs."""
    source = tmp_path / "source.xdmf"
    _write_time_series(source, times)
    requested_time = times[0] + 0.2 * (times[-1] - times[0])
    output = tmp_path / "interpolated.xdmf"
    interpolate_xdmf(source, output, np.array([requested_time]))
    with meshio.xdmf.TimeSeriesReader(output) as reader:
        reader.read_points_cells()
        time, point_data, _ = reader.read_data(0)
    assert time == pytest.approx(requested_time)
    np.testing.assert_allclose(point_data["pressure"], 10 * requested_time)


def test_interpolate_xdmf_preserves_requested_order_and_repeated_times(
    tmp_path: Path, time_series: Path
) -> None:
    """Keep the requested time sequence while validating all values against the range."""
    output = tmp_path / "interpolated.xdmf"
    requested_times = np.array([10.0, 2.0, 2.0, 0.0])
    interpolate_xdmf(time_series, output, requested_times)
    with meshio.xdmf.TimeSeriesReader(output) as reader:
        reader.read_points_cells()
        for step, expected_time in enumerate(requested_times):
            time, point_data, _ = reader.read_data(step)
            assert time == pytest.approx(expected_time)
            np.testing.assert_allclose(point_data["pressure"], 10 * expected_time)


def test_interpolate_xdmf_rejects_in_place_overwrite(time_series: Path) -> None:
    """Do not truncate source data when input and output identify the same file."""
    original = time_series.read_bytes()
    with pytest.raises(ValueError, match="Input and output XDMF paths must differ"):
        interpolate_xdmf(time_series, time_series, np.array([2.0]), overwrite=True)
    assert time_series.read_bytes() == original


def test_vtk_to_xdmf_rejects_empty_directory(tmp_path: Path) -> None:
    """Explain missing VTK input without creating output files."""
    output = tmp_path / "empty.xdmf"
    with pytest.raises(ValueError, match="No VTK files found"):
        vtks_to_xdmf(tmp_path, output)
    assert not output.exists()


def test_vtk_to_xdmf_supports_concurrent_outputs_with_same_filename(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Independent conversions do not share HDF files or change the working directory."""
    monkeypatch.chdir(tmp_path)
    outputs = [
        tmp_path / directory / "simulation.xdmf" for directory in ("first", "second")
    ]
    vtk_dir = RESOURCES_DIR / "vtk" / "vtk_single"
    with ThreadPoolExecutor(max_workers=2) as executor:
        futures = [executor.submit(vtks_to_xdmf, vtk_dir, output) for output in outputs]
        for future in futures:
            future.result()
    assert Path.cwd() == tmp_path
    assert not (tmp_path / "simulation.h5").exists()
    for output in outputs:
        assert XDMFInfo.from_path(output).num_steps == 1
