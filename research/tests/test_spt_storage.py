"""Check that every numerical artifact is written compressed and reads back."""

from pathlib import Path

import numpy as np
import pytest
import xarray as xr

from research.spt.storage import (
    COMPLEVEL,
    LEGACY_NETCDF_ENGINE,
    NETCDF_ENGINE,
    compression_record,
    load_array,
    load_dataset,
    open_dataset,
    save_array,
    save_arrays,
    write_dataset,
)


def _cells(times: int = 200, cells: int = 400) -> xr.Dataset:
    """Build a float32 cell history with the smooth structure of a real run."""
    generator = np.random.default_rng(7)
    axis = np.linspace(0.0, 1.0, cells, dtype=np.float32)
    ramp = np.linspace(0.0, 1.0, times, dtype=np.float32)[:, None]
    smooth = (ramp * axis[None, :]).astype(np.float32)
    return xr.Dataset(
        {
            "concentration": (("time", "cell"), smooth),
            "necrosis": (("time", "cell"), (smooth > 0.5).astype(np.float32)),
            "noise": (
                ("time", "cell"),
                generator.random((times, cells), dtype=np.float32),
            ),
        },
        coords={
            "time": np.arange(times, dtype=float),
            "cell": np.arange(cells, dtype=np.int32),
        },
    )


def test_written_cache_is_compressed_and_round_trips(tmp_path: Path) -> None:
    """Store the arrays deflated and return identical values and dtypes."""
    data = _cells()
    compressed = tmp_path / "cells.nc"
    write_dataset(data, compressed)
    uncompressed = tmp_path / "legacy.nc"
    data.to_netcdf(uncompressed, engine=LEGACY_NETCDF_ENGINE)

    # NetCDF3 cannot compress, which is why the caches moved to this engine.
    assert compressed.stat().st_size < uncompressed.stat().st_size / 2

    result = load_dataset(compressed)
    assert set(result.data_vars) == set(data.data_vars)
    for name, variable in data.data_vars.items():
        assert result[name].dtype == variable.dtype
        np.testing.assert_array_equal(result[name].values, variable.values)
    np.testing.assert_array_equal(result.time.values, data.time.values)


def test_every_array_declares_the_expected_filters(tmp_path: Path) -> None:
    """Apply zlib and the shuffle filter to each stored array, not just the file."""
    h5py = pytest.importorskip("h5py")
    path = tmp_path / "cells.nc"
    write_dataset(_cells(times=20, cells=30), path)
    with h5py.File(path, "r") as handle:
        for name in ("concentration", "necrosis", "noise"):
            dataset = handle[name]
            assert dataset.compression == "gzip"
            assert dataset.compression_opts == COMPLEVEL
            assert dataset.shuffle


def test_scalar_variables_stay_writable(tmp_path: Path) -> None:
    """Skip the filters HDF5 cannot apply to an unchunked scalar."""
    data = xr.Dataset(
        {"total": ((), np.float32(2.5)), "series": (("time",), np.arange(4.0))},
        coords={"time": np.arange(4.0)},
    )
    path = tmp_path / "summary.nc"
    write_dataset(data, path)
    result = load_dataset(path)
    assert float(result["total"]) == pytest.approx(2.5)
    np.testing.assert_array_equal(result["series"].values, np.arange(4.0))


def test_legacy_uncompressed_cache_still_reads(tmp_path: Path) -> None:
    """Keep an artifact written before compression readable until its stage rebuilds."""
    data = _cells(times=4, cells=5)
    path = tmp_path / "legacy.nc"
    data.to_netcdf(path, engine=LEGACY_NETCDF_ENGINE)
    np.testing.assert_array_equal(
        load_dataset(path)["concentration"].values, data["concentration"].values
    )
    with open_dataset(path) as opened:
        np.testing.assert_array_equal(
            opened["concentration"].values, data["concentration"].values
        )


def test_opened_dataset_closes_on_failure(tmp_path: Path) -> None:
    """Release the file handle even when the caller raises inside the context."""
    path = tmp_path / "cells.nc"
    write_dataset(_cells(times=3, cells=4), path)
    with pytest.raises(RuntimeError), open_dataset(path) as data:
        assert data["concentration"].size
        raise RuntimeError("caller failed")
    # A closed handle can be reopened; an escaped one would keep the file locked.
    with open_dataset(path) as reopened:
        assert reopened["concentration"].size


def test_sample_archives_are_compressed_and_round_trip(tmp_path: Path) -> None:
    """Deflate plot samples, which store repetitive face and sector indices."""
    values = np.repeat(np.arange(500, dtype=np.float32), 20)
    sectors = np.repeat(np.arange(6, dtype=np.int32), 1000)
    compressed = tmp_path / "samples.npz"
    save_arrays(compressed, values=values, sectors=sectors)
    plain = tmp_path / "plain.npz"
    np.savez(plain, values=values, sectors=sectors)
    assert compressed.stat().st_size < plain.stat().st_size / 2
    with np.load(compressed, allow_pickle=False) as archive:
        np.testing.assert_array_equal(archive["values"], values)
        np.testing.assert_array_equal(archive["sectors"], sectors)


def test_single_array_round_trips_in_either_layout(tmp_path: Path) -> None:
    """Read a compressed archive and a legacy ``.npy`` through one entry point."""
    heights = np.linspace(0.0, 1.0, 64)
    archive = tmp_path / "heights.npz"
    save_array(archive, heights)
    np.testing.assert_array_equal(load_array(archive), heights)
    legacy = tmp_path / "heights.npy"
    np.save(legacy, heights, allow_pickle=False)
    np.testing.assert_array_equal(load_array(legacy), heights)


def test_reserved_archive_name_is_rejected(tmp_path: Path) -> None:
    """Refuse an array named like the pickle switch instead of silently dropping it."""
    with pytest.raises(ValueError, match="allow_pickle"):
        save_arrays(tmp_path / "samples.npz", allow_pickle=np.arange(3))


def test_compression_record_describes_the_written_settings() -> None:
    """Record settings a reader can compare against a manifest."""
    record = compression_record()
    assert record["netcdf_engine"] == NETCDF_ENGINE
    assert str(COMPLEVEL) in str(record["netcdf_compression"])
    assert str(COMPLEVEL) in str(record["hdf_compression"])
