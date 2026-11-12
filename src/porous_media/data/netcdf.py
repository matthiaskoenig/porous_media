"""Write and read NetCDF datasets with compression decided in one place.

Simulation histories are large: one uncompressed cell dataset of a single run
measured 70 MB where the same float32 arrays with deflate and the byte-shuffle
filter took 13 MB. NetCDF3, which :mod:`scipy` writes, has no compression at
all, so every dataset here goes through the HDF5-backed ``h5netcdf`` engine.

Compression is lossless; the stored dtypes are unchanged. Reading accepts both
layouts so a dataset written before compression stays readable.
"""

from collections.abc import Generator
from contextlib import contextmanager
from pathlib import Path
from typing import Any

import xarray as xr

NETCDF_ENGINE = "h5netcdf"
"""Engine used for writing; NetCDF3 cannot compress."""

LEGACY_NETCDF_ENGINE = "scipy"
"""Engine of uncompressed NetCDF3 files written before compression."""

COMPLEVEL = 4
"""Deflate level. Measured against level 1 on simulation histories, level 4 is
12% smaller for about a third more write time, while levels above 6 gain less
than 0.3%."""


HDF_COMPRESSION: dict[str, Any] = {
    "compression": "gzip",
    "compression_opts": COMPLEVEL,
    "shuffle": True,
}
"""Keyword arguments for :meth:`h5py.Group.create_dataset`, so HDF sidecars use
the same filters as the NetCDF datasets."""


def netcdf_encoding(dataset: xr.Dataset) -> dict[str, dict[str, Any]]:
    """Request deflate and shuffle for every array HDF5 can chunk.

    HDF5 applies a filter only to a chunked dataset, so a zero-dimensional or
    empty variable is left out: asking for compression on one is an error rather
    than a no-op. Non-numeric variables are also left out.
    """
    encoding: dict[str, dict[str, Any]] = {}
    for name, variable in dataset.variables.items():
        if variable.ndim == 0 or variable.size == 0:
            continue
        if variable.dtype.kind not in "biufc":
            continue
        encoding[str(name)] = {
            "zlib": True,
            "complevel": COMPLEVEL,
            "shuffle": True,
        }
    return encoding


def write_netcdf(dataset: xr.Dataset, path: Path) -> None:
    """Write one dataset compressed, leaving the caller's object unchanged."""
    dataset.to_netcdf(path, engine=NETCDF_ENGINE, encoding=netcdf_encoding(dataset))


def load_netcdf(path: Path, **kwargs: Any) -> xr.Dataset:
    """Read a dataset into memory, accepting a legacy uncompressed NetCDF3 file."""
    try:
        return xr.load_dataset(path, engine=NETCDF_ENGINE, **kwargs)
    except (OSError, ValueError):
        return xr.load_dataset(path, engine=LEGACY_NETCDF_ENGINE, **kwargs)


@contextmanager
def open_netcdf(path: Path, **kwargs: Any) -> Generator[xr.Dataset]:
    """Open a dataset lazily, accepting a legacy uncompressed NetCDF3 file."""
    try:
        data = xr.open_dataset(path, engine=NETCDF_ENGINE, **kwargs)
    except (OSError, ValueError):
        data = xr.open_dataset(path, engine=LEGACY_NETCDF_ENGINE, **kwargs)
    try:
        yield data
    finally:
        data.close()
