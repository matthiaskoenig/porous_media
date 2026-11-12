"""Test xdmf functionality."""

from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import h5py
import meshio
import numpy as np
import pytest
from meshio._exceptions import WriteError

from porous_media import RESOURCES_DIR
from porous_media.data.xdmf_tools import (
    AttributeType,
    DataLimits,
    XDMFInfo,
    XDMFTimeSeriesReader,
    XDMFTimeSeriesWriter,
    interpolate_xdmf,
    vtks_to_xdmf,
    xdmfs_from_directory,
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


def test_xdmf_writer_compresses_and_shares_repeated_arrays(tmp_path: Path) -> None:
    """Deflate each array and store a repeated time-invariant field only once."""
    output = tmp_path / "simulation.xdmf"
    points = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    fixed = np.linspace(0.0, 1.0, 3)
    with XDMFTimeSeriesWriter(output, share_identical_arrays=True) as writer:
        writer.write_points_cells(points, [("triangle", np.array([[0, 1, 2]]))])
        for step in range(4):
            writer.write_data(
                float(step),
                # "zonation" never changes, so every step repeats one array.
                point_data={"zonation": fixed, "pressure": fixed + step},
            )

    with h5py.File(output.with_suffix(".h5"), "r") as handle:
        names = list(handle)
        assert names
        for name in names:
            assert handle[name].compression == "gzip"
            assert handle[name].shuffle
        # Identical arrays share one HDF object instead of being stored again.
        identities = {handle[name].id.__hash__() for name in names}
        assert len(identities) < len(names)

    with meshio.xdmf.TimeSeriesReader(output) as reader:
        reader.read_points_cells()
        for step in range(4):
            _, point_data, _ = reader.read_data(step)
            np.testing.assert_allclose(point_data["zonation"], fixed)
            np.testing.assert_allclose(point_data["pressure"], fixed + step)


def test_xdmf_writer_keeps_arrays_independent_by_default(tmp_path: Path) -> None:
    """Store each array separately so a sidecar stays safe to edit in place."""
    output = tmp_path / "independent.xdmf"
    fixed = np.linspace(0.0, 1.0, 3)
    with XDMFTimeSeriesWriter(output) as writer:
        writer.write_points_cells(
            np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
            [("triangle", np.array([[0, 1, 2]]))],
        )
        for step in range(3):
            writer.write_data(float(step), point_data={"zonation": fixed})

    with h5py.File(output.with_suffix(".h5"), "r") as handle:
        names = list(handle)
        assert len({handle[name].id.__hash__() for name in names}) == len(names)
        for name in names:
            assert handle[name].compression == "gzip"

    # Editing one step in place must leave the other steps alone.
    with meshio.xdmf.TimeSeriesReader(output) as reader:
        reader.read_points_cells()
        target = reader.read_data(1)
    assert target is not None
    with h5py.File(output.with_suffix(".h5"), "r+") as handle:
        edited = [
            name
            for name in handle
            if handle[name].shape == fixed.shape
            and np.allclose(handle[name][...], fixed)
        ][1]
        handle[edited][...] = np.full(3, 9.0)
    with meshio.xdmf.TimeSeriesReader(output) as reader:
        reader.read_points_cells()
        changed = sum(
            bool(np.allclose(reader.read_data(step)[1]["zonation"], 9.0))
            for step in range(3)
        )
    assert changed == 1


def test_xdmf_writer_keeps_xml_data_format_uncompressed(tmp_path: Path) -> None:
    """Leave inline XML output alone; only the HDF sidecar carries filters."""
    output = tmp_path / "inline.xdmf"
    with XDMFTimeSeriesWriter(output, data_format="XML") as writer:
        writer.write_points_cells(
            np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
            [("triangle", np.array([[0, 1, 2]]))],
        )
        writer.write_data(0.0, point_data={"pressure": np.array([1.0, 2.0, 3.0])})
    assert not output.with_suffix(".h5").exists()
    with meshio.xdmf.TimeSeriesReader(output) as reader:
        reader.read_points_cells()
        _, point_data, _ = reader.read_data(0)
    np.testing.assert_array_equal(point_data["pressure"], [1.0, 2.0, 3.0])


def _write_series_steps(path: Path, steps: int, *, stacked: bool) -> None:
    """Write a fixed-mesh series whose field changes only slightly per step."""
    count = 64
    points = np.zeros((count + 2, 3))
    points[:, 0] = np.arange(count + 2)
    points[1::2, 1] = 1.0
    triangles = np.array([[i, i + 1, i + 2] for i in range(count)])
    base = np.linspace(0.0, 1.0, count, dtype=np.float32)
    with XDMFTimeSeriesWriter(path, stack_time_steps=stacked) as writer:
        writer.write_points_cells(points, [("triangle", triangles)])
        for step in range(steps):
            writer.write_data(
                float(step),
                point_data={"slow": np.linspace(0.0, 1.0, count + 2) + step * 1e-4},
                cell_data={"field": [base + step]},
            )


def test_stacked_series_round_trips_through_the_matching_reader(
    tmp_path: Path,
) -> None:
    """Resolve each step's HyperSlab back to the values that were written."""
    path = tmp_path / "stacked.xdmf"
    _write_series_steps(path, 12, stacked=True)
    base = np.linspace(0.0, 1.0, 64, dtype=np.float32)

    with h5py.File(path.with_suffix(".h5"), "r") as handle:
        # One growing dataset per attribute, chunked across several steps.
        assert handle["Cell/field"].shape == (12, 64)
        assert handle["Node/slow"].shape == (12, 66)
        assert handle["Cell/field"].chunks[0] > 1
        assert handle["Cell/field"].compression == "gzip"

    with XDMFTimeSeriesReader(path) as reader:
        reader.read_points_cells()
        assert reader.num_steps == 12
        for step in range(12):
            time, point_data, cell_data = reader.read_data(step)
            assert time == pytest.approx(step)
            np.testing.assert_allclose(np.concatenate(cell_data["field"]), base + step)
            np.testing.assert_allclose(
                point_data["slow"], np.linspace(0.0, 1.0, 66) + step * 1e-4
            )


def test_stacking_shrinks_a_slowly_changing_series(tmp_path: Path) -> None:
    """Give the compressor consecutive steps instead of one frame at a time."""
    per_step = tmp_path / "per-step.xdmf"
    stacked = tmp_path / "stacked.xdmf"
    _write_series_steps(per_step, 300, stacked=False)
    _write_series_steps(stacked, 300, stacked=True)
    assert (
        stacked.with_suffix(".h5").stat().st_size
        < per_step.with_suffix(".h5").stat().st_size / 2
    )


def test_matching_reader_also_reads_an_unstacked_series(tmp_path: Path) -> None:
    """Stay a drop-in replacement for the plain reader."""
    path = tmp_path / "per-step.xdmf"
    _write_series_steps(path, 4, stacked=False)
    with XDMFTimeSeriesReader(path) as ours:
        ours.read_points_cells()
        with meshio.xdmf.TimeSeriesReader(path) as theirs:
            theirs.read_points_cells()
            for step in range(4):
                _, _, mine = ours.read_data(step)
                _, _, reference = theirs.read_data(step)
                np.testing.assert_array_equal(
                    np.concatenate(mine["field"]),
                    np.concatenate(reference["field"]),
                )


def test_stacked_series_rejects_a_changing_field_shape(tmp_path: Path) -> None:
    """Refuse a series whose attribute changes shape rather than writing a wrong slab."""
    path = tmp_path / "changing.xdmf"
    points = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    with (
        pytest.raises(WriteError, match="shape or type"),
        XDMFTimeSeriesWriter(path, stack_time_steps=True) as writer,
    ):
        writer.write_points_cells(points, [("triangle", np.array([[0, 1, 2]]))])
        writer.write_data(0.0, point_data={"p": np.zeros(3)})
        writer.write_data(1.0, point_data={"p": np.zeros((3, 3))})


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


def test_directory_conversion_rejects_flattened_name_collisions(tmp_path: Path) -> None:
    """Detect ambiguous output names before converting either input directory."""
    inputs = tmp_path / "inputs"
    for relative in ("a/b", "a__b"):
        directory = inputs / relative
        directory.mkdir(parents=True)
        (directory / "time.vtk").write_text("unused", encoding="utf-8")
    output = tmp_path / "output"
    with pytest.raises(ValueError, match="same output"):
        xdmfs_from_directory(inputs, output)
    assert not output.exists()


def test_conversion_preserves_orphaned_sidecar(tmp_path: Path) -> None:
    """An unrelated HDF file must not be replaced when only the XDMF file is absent."""
    output = tmp_path / "source.xdmf"
    sidecar = output.with_suffix(".h5")
    sidecar.write_bytes(b"unrelated data")
    with pytest.raises(FileExistsError, match="sidecar"):
        vtks_to_xdmf(RESOURCES_DIR / "vtk" / "vtk_single", output)
    assert sidecar.read_bytes() == b"unrelated data"
    assert not output.exists()


def test_limits_merge_supports_different_field_sets() -> None:
    """Merge heterogeneous simulations without requiring the first to contain every field."""
    first = DataLimits({"pressure": (1.0, 2.0)})
    result = DataLimits.merge_limits(
        [first, DataLimits({"pressure": (0.0, 3.0), "new": (4.0, 5.0)})]
    )
    assert result.limits == {"pressure": (0.0, 3.0), "new": (4.0, 5.0)}
    assert first.limits == {"pressure": (1.0, 2.0)}


def test_conversion_refreshes_an_existing_limits_cache(tmp_path: Path) -> None:
    """Recalculate limits whenever a new time series is generated."""
    output = tmp_path / "simulation.xdmf"
    cache = DataLimits.json_path_from_xdmf(output)
    cache.write_text(
        '{"limits": {"effective_fluid_pressure_TPM": [-999, -999]}}', encoding="utf-8"
    )
    vtks_to_xdmf(RESOURCES_DIR / "vtk" / "vtk_single", output)
    limits = DataLimits.from_xdmf(output)
    assert limits.limits["effective_fluid_pressure_TPM"] != (-999, -999)


def test_conversion_rejects_symlinked_sidecar(tmp_path: Path) -> None:
    """Preserve a symlink target even when overwrite is explicitly enabled."""
    output = tmp_path / "simulation.xdmf"
    target = tmp_path / "unrelated.h5"
    target.write_bytes(b"unrelated data")
    try:
        output.with_suffix(".h5").symlink_to(target)
    except OSError:
        pytest.skip("Creating symbolic links is unavailable")
    with pytest.raises(ValueError, match="symbolic links"):
        vtks_to_xdmf(RESOURCES_DIR / "vtk" / "vtk_single", output, overwrite=True)
    assert target.read_bytes() == b"unrelated data"


def test_xdmf_writer_closes_hdf_after_metadata_failure(tmp_path: Path) -> None:
    """Allow reopening a sidecar after XML serialization fails."""
    output = tmp_path / "simulation.xdmf"
    output.mkdir()
    writer = XDMFTimeSeriesWriter(output)
    with pytest.raises(OSError), writer:
        writer.write_points_cells(np.eye(3), [("triangle", np.array([[0, 1, 2]]))])
    with h5py.File(output.with_suffix(".h5"), "w") as reopened:
        reopened.create_dataset("valid", data=[1])


def test_xdmf_metadata_preserves_six_component_tensor_shapes(tmp_path: Path) -> None:
    """Describe symmetric tensors on points and cells without artificial dimensions."""
    source = tmp_path / "tensor.xdmf"
    with meshio.xdmf.TimeSeriesWriter(source, data_format="XML") as writer:
        writer.write_points_cells(np.eye(3), [("triangle", np.array([[0, 1, 2]]))])
        writer.write_data(
            0.0,
            point_data={"point_stress": np.ones((3, 6))},
            cell_data={"stress": [np.ones((1, 6))]},
        )
    info = XDMFInfo.from_path(source)
    assert info.point_data["point_stress"].attribute_type == AttributeType.TENSOR6
    assert info.cell_data["stress"].attribute_type == AttributeType.TENSOR6
    assert info.point_data["point_stress"].shape == (3, 6)
    assert info.cell_data["stress"].shape == (1, 6)
