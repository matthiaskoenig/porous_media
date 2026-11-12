"""Check lossless NetCDF compression and the HDF settings shared with XDMF."""

from pathlib import Path

import h5py
import numpy as np
import pytest
import xarray as xr

from porous_media.data.netcdf import (
    COMPLEVEL,
    HDF_COMPRESSION,
    LEGACY_NETCDF_ENGINE,
    load_netcdf,
    netcdf_encoding,
    open_netcdf,
    write_netcdf,
)


def _history(times: int = 100, cells: int = 300) -> xr.Dataset:
    """Build a float32 history with the smooth structure of a simulated field."""
    axis = np.linspace(0.0, 1.0, cells, dtype=np.float32)
    ramp = np.linspace(0.0, 1.0, times, dtype=np.float32)[:, None]
    values = (ramp * axis[None, :]).astype(np.float32)
    return xr.Dataset(
        {
            "pressure": (("time", "cell"), values),
            "necrosis": (("time", "cell"), (values > 0.5).astype(np.float32)),
        },
        coords={
            "time": np.arange(times, dtype=float),
            "cell": np.arange(cells, dtype=np.int32),
        },
    )


def test_compression_is_lossless_and_much_smaller(tmp_path: Path) -> None:
    """Return identical values and dtypes from a far smaller file."""
    data = _history()
    compressed = tmp_path / "history.nc"
    uncompressed = tmp_path / "history-netcdf3.nc"
    write_netcdf(data, compressed)
    data.to_netcdf(uncompressed, engine=LEGACY_NETCDF_ENGINE)
    assert compressed.stat().st_size < uncompressed.stat().st_size / 2

    result = load_netcdf(compressed)
    for name, variable in data.data_vars.items():
        assert result[name].dtype == variable.dtype
        np.testing.assert_array_equal(result[name].values, variable.values)


def test_encoding_covers_arrays_and_skips_unchunkable_variables() -> None:
    """Request filters for each array, and omit the scalars HDF5 cannot chunk."""
    data = xr.Dataset(
        {
            "series": (("time",), np.arange(4.0)),
            "total": ((), np.float32(1.5)),
            "empty": (("missing",), np.empty(0)),
            "label": (("tag",), np.array(["a", "b"])),
        },
        coords={"time": np.arange(4.0), "tag": np.array(["x", "y"])},
    )
    encoding = netcdf_encoding(data)
    assert encoding["series"] == {
        "zlib": True,
        "complevel": COMPLEVEL,
        "shuffle": True,
    }
    for skipped in ("total", "empty", "label"):
        assert skipped not in encoding


def test_unchunkable_variables_still_write(tmp_path: Path) -> None:
    """Write a dataset holding a scalar without asking HDF5 for an invalid filter."""
    data = xr.Dataset(
        {"total": ((), np.float32(2.5)), "series": (("time",), np.arange(3.0))},
        coords={"time": np.arange(3.0)},
    )
    path = tmp_path / "summary.nc"
    write_netcdf(data, path)
    result = load_netcdf(path)
    assert float(result["total"]) == pytest.approx(2.5)


def test_legacy_netcdf3_file_still_reads(tmp_path: Path) -> None:
    """Read a file written before compression, loaded and opened."""
    data = _history(times=3, cells=4)
    path = tmp_path / "legacy.nc"
    data.to_netcdf(path, engine=LEGACY_NETCDF_ENGINE)
    np.testing.assert_array_equal(
        load_netcdf(path)["pressure"].values, data["pressure"].values
    )
    with open_netcdf(path) as opened:
        np.testing.assert_array_equal(
            opened["pressure"].values, data["pressure"].values
        )


def test_opened_dataset_closes_when_the_caller_raises(tmp_path: Path) -> None:
    """Release the handle on failure so the file can be reopened."""
    path = tmp_path / "history.nc"
    write_netcdf(_history(times=3, cells=4), path)
    with pytest.raises(RuntimeError), open_netcdf(path) as data:
        assert data["pressure"].size
        raise RuntimeError("caller failed")
    with open_netcdf(path) as reopened:
        assert reopened["pressure"].size


def test_hdf_settings_are_accepted_by_h5py(tmp_path: Path) -> None:
    """Keep the shared settings usable as `create_dataset` keywords."""
    path = tmp_path / "arrays.h5"
    values = np.arange(1000, dtype=np.float32)
    with h5py.File(path, "w") as handle:
        handle.create_dataset("values", data=values, **HDF_COMPRESSION)
    with h5py.File(path, "r") as handle:
        dataset = handle["values"]
        assert dataset.compression == "gzip"
        assert dataset.compression_opts == COMPLEVEL
        assert dataset.shuffle
        np.testing.assert_array_equal(dataset[...], values)
