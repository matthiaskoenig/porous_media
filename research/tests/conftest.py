"""Shared guards for the research test suite."""

from collections.abc import Callable
from pathlib import Path

import h5py
import pytest

from porous_media.data.xdmf_tools import XDMFTimeSeriesReader


@pytest.fixture
def reject_history_hdf(monkeypatch: pytest.MonkeyPatch) -> Callable[[str], None]:
    """Fail a test that opens a raw or reconstructed XDMF history.

    Several stages must work from the compact local artifacts rather than
    reading a full history again. Those histories are the HDF sidecars of an
    XDMF time series and always carry the ``.h5`` suffix, while the caches and
    published results are NetCDF (``.nc``). Banning :class:`h5py.File` outright
    used to express that, because an uncompressed NetCDF3 file was read through
    :mod:`scipy` rather than HDF5; now that those datasets are written
    compressed they are HDF5 too, so the suffix is what separates them.

    The time series reader is barred as well, on the class rather than on the
    :mod:`meshio` name it inherits from, so that patching holds however the
    reader was imported.
    """

    def install(message: str) -> None:
        original = h5py.File

        def guard(name: object, *args: object, **kwargs: object) -> object:
            if Path(str(name)).suffix == ".nc":
                return original(name, *args, **kwargs)
            pytest.fail(f"{message}: {name}")

        def reader(*args: object, **kwargs: object) -> None:
            pytest.fail(message)

        monkeypatch.setattr("h5py.File", guard)
        monkeypatch.setattr(XDMFTimeSeriesReader, "__init__", reader)

    return install
