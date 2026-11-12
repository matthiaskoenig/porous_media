"""Write every numerical SPT artifact compressed, and read it back by name.

The SPT pipeline stores large per-cell histories, and written without
compression these dominate the results tree: one uncompressed cell cache
measured 419 MB against 36 MB for the same float32 arrays deflated. Every cache
therefore goes through this module, so the engine, the compression level and the
read path are decided in one place.

NetCDF handling lives in :mod:`porous_media.data.netcdf` because the library
writes the same kind of dataset; this module adds the archive helpers for the
sampled plot data and the manifest record. Compression is lossless and leaves
the stored dtypes unchanged. Reading accepts the uncompressed layout too, so an
artifact written before compression stays readable until the stage that owns it
rebuilds.
"""

from collections.abc import Mapping
from pathlib import Path
from typing import Any, cast

import numpy as np

from porous_media.data.netcdf import COMPLEVEL as COMPLEVEL
from porous_media.data.netcdf import HDF_COMPRESSION as HDF_COMPRESSION
from porous_media.data.netcdf import LEGACY_NETCDF_ENGINE as LEGACY_NETCDF_ENGINE
from porous_media.data.netcdf import NETCDF_ENGINE as NETCDF_ENGINE
from porous_media.data.netcdf import load_netcdf as load_dataset
from porous_media.data.netcdf import open_netcdf as open_dataset
from porous_media.data.netcdf import write_netcdf as write_dataset

__all__ = [
    "COMPLEVEL",
    "HDF_COMPRESSION",
    "LEGACY_NETCDF_ENGINE",
    "NETCDF_ENGINE",
    "compression_record",
    "load_array",
    "load_dataset",
    "open_dataset",
    "save_array",
    "save_arrays",
    "write_dataset",
]


def save_arrays(path: Path, **arrays: np.ndarray) -> None:
    """Write a compressed ``.npz`` archive of plot samples."""
    if "allow_pickle" in arrays:
        raise ValueError("'allow_pickle' is reserved and cannot name an array")
    np.savez_compressed(path, allow_pickle=False, **cast(dict[str, Any], arrays))


def save_array(path: Path, array: np.ndarray) -> None:
    """Write one array as a compressed single-entry ``.npz`` archive."""
    save_arrays(path, array=array)


def load_array(path: Path) -> np.ndarray:
    """Read an array written by :func:`save_array`, accepting a legacy ``.npy``."""
    if path.suffix == ".npy":
        return np.load(path, allow_pickle=False)
    with np.load(path, allow_pickle=False) as archive:
        return archive["array"]


def compression_record() -> Mapping[str, Any]:
    """Describe the storage settings for an artifact manifest."""
    return {
        "netcdf_engine": NETCDF_ENGINE,
        "netcdf_compression": f"zlib level {COMPLEVEL} with shuffle",
        "array_compression": "deflate (savez_compressed)",
        "hdf_compression": f"gzip level {COMPLEVEL} with shuffle",
    }
