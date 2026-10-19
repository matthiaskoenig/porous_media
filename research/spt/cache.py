"""Build local SPT caches with provenance, bounded reads, and atomic publication."""

import json
import os
import shutil
import tempfile
import time
import uuid
from collections.abc import Generator
from contextlib import ExitStack, contextmanager
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Literal, cast

import h5py
import numpy as np
import xarray as xr

from research.spt.catalog import (
    CatalogValue,
    ParsedRun,
    RunSource,
    hdf_reference,
    inspect_run,
)

type Profile = Literal["analysis", "full"]
type LegacyPolicy = Literal["sample", "full", "ignore"]
CACHE_VERSION = 1
ANALYSIS_FIELDS = (
    "rr_(S_ext)",
    "rr_(P_ext)",
    "rr_protein",
    "rr_(T)",
    "rr_necrosis",
    "element_volume_point_TPM",
    "rr_position",
    "rr_Vext",
    "rr_Vli",
    "pressure",
    "rr_zonation_pattern",
    "cell_type",
)
SUMMARY_FIELDS = ANALYSIS_FIELDS[:5]


@dataclass
class ProcessingResult:
    """Report per-run success, failure, and cache reuse to the coordinator."""

    study: str
    run_id: str
    status: str
    message: str
    elapsed_s: float
    row: dict[str, CatalogValue]


def validate_roots(input_root: Path, output_root: Path) -> tuple[Path, Path]:
    """Keep all derived outputs separate from the input tree."""
    source = input_root.expanduser().resolve()
    target = output_root.expanduser().resolve()
    if target.is_relative_to(source) or source.is_relative_to(target):
        raise ValueError(
            "Input and output directories must be separate, nonoverlapping trees"
        )
    return source, target


def run_directory(source: RunSource, output_root: Path) -> Path:
    """Locate the local cache without depending on input file ordering."""
    path = output_root.expanduser().resolve() / source.study / "runs" / source.run_id
    if path.resolve() != path.absolute():
        raise ValueError("Cache directories must not traverse symbolic links")
    return path


def fingerprint(path: Path) -> dict[str, int | str | None]:
    """Record cheap change detection without hashing whole network files."""
    try:
        stat = path.stat()
    except FileNotFoundError:
        return {"path": str(path), "size": None, "mtime_ns": None}
    return {"path": str(path), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns}


def _configuration(
    source: RunSource, profile: Profile, legacy: LegacyPolicy
) -> dict[str, object]:
    return {
        "version": CACHE_VERSION,
        "profile": profile,
        "legacy_policy": legacy,
        "analysis_fields": list(ANALYSIS_FIELDS) if profile == "analysis" else None,
        "condition": asdict(source.condition) if source.condition is not None else None,
    }


def _manifest(source: RunSource, output_root: Path) -> dict[str, object] | None:
    path = run_directory(source, output_root) / "manifest.json"
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (FileNotFoundError, ValueError):
        return None
    return value if isinstance(value, dict) else None


def _generation(directory: Path, manifest: dict[str, object]) -> Path:
    name = manifest.get("generation")
    if not isinstance(name, str) or not name.startswith("generation-"):
        raise ValueError("Invalid cache generation")
    uuid.UUID(name.removeprefix("generation-"))
    path = directory / name
    if path.is_symlink() or not path.is_dir():
        raise ValueError("Cache generation missing or a symbolic link")
    return path


def _unchanged(records: object, root: Path) -> bool:
    if not isinstance(records, list) or not records:
        return False
    for record in records:
        if not isinstance(record, dict) or not isinstance(record.get("path"), str):
            return False
        path = Path(record["path"])
        if not path.resolve().is_relative_to(root) or fingerprint(path) != record:
            return False
    return True


def cache_status(
    source: RunSource, output_root: Path, profile: Profile, legacy: LegacyPolicy
) -> tuple[str, dict[str, object] | None]:
    """Check a completed cache using file metadata, without reading arrays or XML."""
    manifest = _manifest(source, output_root)
    if manifest is None:
        return "missing", None
    try:
        directory = run_directory(source, output_root)
        generation = _generation(directory, manifest)
        if manifest.get("configuration") != _configuration(source, profile, legacy):
            return "configuration changed", manifest
        if not _unchanged(manifest.get("sources"), source.root):
            return "source changed", manifest
        if not _unchanged(manifest.get("outputs"), generation):
            return "output changed", manifest
        if not isinstance(manifest.get("catalog"), dict):
            return "invalid catalog", manifest
    except (OSError, ValueError):
        return "invalid", manifest
    return "current", manifest


def cached_catalog(
    source: RunSource, output_root: Path
) -> dict[str, CatalogValue] | None:
    """Reuse catalogue metadata independently of the numerical cache profile."""
    manifest = _manifest(source, output_root)
    if manifest is None or not _unchanged(manifest.get("sources"), source.root):
        return None
    configuration = manifest.get("configuration")
    condition = asdict(source.condition) if source.condition is not None else None
    if (
        not isinstance(configuration, dict)
        or configuration.get("version") != CACHE_VERSION
        or configuration.get("condition") != condition
    ):
        return None
    row = manifest.get("catalog")
    return cast(dict[str, CatalogValue], row) if isinstance(row, dict) else None


@contextmanager
def _run_lock(directory: Path) -> Generator[None]:
    """Use an OS lock that is released automatically when a worker exits."""
    directory.mkdir(parents=True, exist_ok=True)
    if directory.is_symlink() or (directory / ".lock").is_symlink():
        raise ValueError("Cache directories and locks must not be symbolic links")
    with (directory / ".lock").open("a+b") as handle:
        if os.name == "nt":
            import msvcrt

            handle.write(b"\0")
            handle.flush()
            handle.seek(0)
            msvcrt.locking(handle.fileno(), msvcrt.LK_NBLCK, 1)
        else:
            import fcntl

            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        try:
            yield
        finally:
            if os.name == "nt":
                handle.seek(0)
                msvcrt.locking(handle.fileno(), msvcrt.LK_UNLCK, 1)
            else:
                fcntl.flock(handle, fcntl.LOCK_UN)


def _fields(parsed: ParsedRun, profile: Profile) -> dict[str, tuple[str, ...]]:
    if profile == "full":
        return parsed.fields
    missing = set(ANALYSIS_FIELDS) - set(parsed.fields["Cell"])
    if missing:
        raise ValueError(f"Missing analysis fields: {', '.join(sorted(missing))}")
    return {"Cell": ANALYSIS_FIELDS, "Node": ()}


def _read_frame(
    parsed: ParsedRun,
    step: int,
    center: str,
    names: tuple[str, ...],
    handles: dict[str, h5py.File],
) -> dict[str, np.ndarray]:
    result = {}
    attrs = {
        attr.attrib["Name"]: attr
        for attr in parsed.grids[step].findall("Attribute")
        if attr.get("Center", "Node") == center
    }
    for name in names:
        item = attrs[name].find("DataItem")
        if item is None:
            raise ValueError(f"Missing DataItem for {name}")
        filename, reference = hdf_reference(item)
        dataset = handles[filename][reference]
        if not isinstance(dataset, h5py.Dataset):
            raise ValueError(f"Not an HDF dataset: {reference}")
        values = np.asarray(dataset[()])
        shape = tuple(int(value) for value in item.attrib["Dimensions"].split())
        if values.shape != shape:
            raise ValueError(f"HDF shape differs from descriptor for {name}")
        result[name] = values.reshape(-1)
    return result


def _from_hdf(
    parsed: ParsedRun, fields: dict[str, tuple[str, ...]], handles: dict[str, h5py.File]
) -> dict[str, xr.Dataset]:
    datasets = {}
    for center, names in fields.items():
        dimension = "cell" if center == "Cell" else "point"
        count = parsed.count(center)
        arrays: dict[str, np.ndarray] = {}
        for step in range(len(parsed.times)):
            for name, values in _read_frame(
                parsed, step, center, names, handles
            ).items():
                if step == 0:
                    arrays[name] = np.empty(
                        (len(parsed.times), count), dtype=values.dtype
                    )
                elif not np.can_cast(values.dtype, arrays[name].dtype, casting="safe"):
                    arrays[name] = arrays[name].astype(
                        np.result_type(arrays[name].dtype, values.dtype)
                    )
                arrays[name][step] = values
        datasets[center] = xr.Dataset(
            {name: (("time", dimension), values) for name, values in arrays.items()},
            coords={"time": parsed.times, dimension: np.arange(count)},
        )
    return datasets


def _from_legacy(
    parsed: ParsedRun,
    fields: dict[str, tuple[str, ...]],
    handles: dict[str, h5py.File],
    policy: LegacyPolicy,
) -> dict[str, xr.Dataset]:
    datasets = {}
    for center, names in fields.items():
        dimension = "cell" if center == "Cell" else "point"
        count = parsed.count(center)
        if not names:
            datasets[center] = xr.Dataset(
                coords={"time": parsed.times, dimension: np.arange(count)}
            )
            continue
        path = parsed.source.xdmf.with_name(f"{parsed.source.run_id}_{dimension}s.nc")
        with xr.open_dataset(path, engine="scipy", decode_times=False) as data:
            if data.sizes.get(dimension) != count or not np.array_equal(
                data.time.values, parsed.times
            ):
                raise ValueError("Legacy cache dimensions or times differ from XDMF")
            if any(data[name].dims != ("time", dimension) for name in names):
                raise ValueError("Legacy cache field dimensions differ from XDMF")
            selected = (
                data[list(names)].assign_coords({dimension: np.arange(count)}).load()
            )
            datasets[center] = selected.copy(deep=True)
        steps = (
            range(len(parsed.times))
            if policy == "full"
            else sorted({0, len(parsed.times) // 2, len(parsed.times) - 1})
        )
        for step in steps:
            for name, values in _read_frame(
                parsed, step, center, names, handles
            ).items():
                if not np.array_equal(
                    values, datasets[center][name].values[step], equal_nan=True
                ):
                    raise ValueError(
                        f"Legacy cache values differ for {name} at step {step}"
                    )
    return datasets


def summarize(cells: xr.Dataset) -> xr.Dataset:
    """Compute native-time volume-weighted summaries with float64 accumulators.

    Interpolation of these summaries is not equivalent to interpolating changing
    cell values and volumes first. Retain per-cell data for that workflow.
    """
    fields = [name for name in SUMMARY_FIELDS if name in cells]
    if "element_volume_point_TPM" not in cells or not fields:
        return xr.Dataset(coords={"time": cells.time})
    weights = cells["element_volume_point_TPM"].astype("float64")
    if (
        not np.isfinite(weights.values).all()
        or bool((weights < 0).any())
        or bool((weights.sum("cell") <= 0).any())
    ):
        raise ValueError(
            "Cell volumes must be finite, nonnegative, and have positive totals"
        )
    values = {}
    for name in fields:
        data = cells[name].astype("float64")
        if not np.isfinite(data.values).all():
            raise ValueError(f"Nonfinite analysis values: {name}")
        weighted = data.weighted(weights)
        values[f"{name}_mean"] = weighted.mean("cell")
        values[f"{name}_std"] = weighted.std("cell")
    result = xr.Dataset(values)
    result.attrs["definition"] = (
        "Native-time volume-weighted population mean/std; no temporal interpolation"
    )
    return result


def _build(
    source: RunSource,
    output_root: Path,
    profile: Profile,
    legacy: LegacyPolicy,
    force: bool,
) -> tuple[str, dict[str, CatalogValue]]:
    directory = run_directory(source, output_root)
    with _run_lock(directory):
        status, previous = cache_status(source, output_root, profile, legacy)
        if status == "current" and not force and previous is not None:
            return "cached", cast(dict[str, CatalogValue], previous["catalog"])
        parsed = inspect_run(source)
        fields = _fields(parsed, profile)
        paths = list(parsed.sources)
        legacy_paths = [
            source.xdmf.with_name(f"{source.run_id}_{kind}s.nc")
            for kind in ("cell", "point")
        ]
        paths.extend(legacy_paths)
        before = [fingerprint(path) for path in paths]
        method = "HDF source"
        fallback = ""
        with ExitStack() as stack:
            handles = {
                name: stack.enter_context(h5py.File(path, "r"))
                for name, path in parsed.references.items()
            }
            if (
                legacy != "ignore"
                and not force
                and previous is None
                and source.study.endswith("_sixth")
            ):
                try:
                    datasets = _from_legacy(parsed, fields, handles, legacy)
                    method = f"legacy cache ({legacy} validation)"
                except (OSError, ValueError, KeyError, TypeError) as error:
                    fallback = f"; legacy rejected: {error}"
                    datasets = _from_hdf(parsed, fields, handles)
            else:
                datasets = _from_hdf(parsed, fields, handles)
        summaries = summarize(datasets["Cell"])
        if before != [fingerprint(path) for path in paths]:
            raise RuntimeError(
                "Inputs changed during processing; retry when upload is complete"
            )
        with tempfile.TemporaryDirectory(
            prefix=".partial-", dir=directory
        ) as temporary:
            staging = Path(temporary)
            for center, filename in (("Cell", "cells.nc"), ("Node", "points.nc")):
                datasets[center].to_netcdf(staging / filename, engine="scipy")
            summaries.to_netcdf(staging / "summary.nc", engine="scipy")
            limits = {}
            for center, data in datasets.items():
                limits[center] = {}
                for name, variable in data.data_vars.items():
                    values = variable.values
                    finite = values[np.isfinite(values)]
                    limits[center][name] = {
                        "min": float(finite.min()) if finite.size else None,
                        "max": float(finite.max()) if finite.size else None,
                        "nonfinite": int(values.size - finite.size),
                    }
            (staging / "limits.json").write_text(
                json.dumps(limits, indent=2, allow_nan=False) + "\n", encoding="utf-8"
            )
            generation = directory / f"generation-{uuid.uuid4()}"
            staging.rename(generation)
            try:
                if before != [fingerprint(path) for path in paths]:
                    raise RuntimeError(
                        "Inputs changed before publication; retry after upload completes"
                    )
                row = {**parsed.row, "cache_method": method}
                manifest = {
                    "configuration": _configuration(source, profile, legacy),
                    "generation": generation.name,
                    "sources": before,
                    "outputs": [
                        fingerprint(path) for path in sorted(generation.iterdir())
                    ],
                    "catalog": row,
                }
                pointer = directory / f".manifest-{uuid.uuid4()}.json"
                try:
                    pointer.write_text(
                        json.dumps(manifest, indent=2, allow_nan=False) + "\n",
                        encoding="utf-8",
                    )
                    pointer.replace(directory / "manifest.json")
                finally:
                    pointer.unlink(missing_ok=True)
            except BaseException:
                shutil.rmtree(generation)
                raise
        if previous is not None:
            try:
                old = _generation(directory, previous)
                if old != generation:
                    shutil.rmtree(old)
            except (OSError, ValueError):
                pass  # Publication succeeded; an older generation can be cleaned later.
        return method + fallback, row


def process_run(
    source: RunSource,
    output_root: Path,
    profile: Profile = "analysis",
    legacy: LegacyPolicy = "sample",
    force: bool = False,
) -> ProcessingResult:
    """Process one run without allowing one failure to discard other completed runs."""
    start = time.perf_counter()
    try:
        if source.study.endswith("_lobulus"):
            output_root = output_root.expanduser().resolve()
            directory = run_directory(source, output_root)
            if (
                source.root.resolve() != output_root
                or not source.xdmf.resolve().is_relative_to(
                    output_root / source.study / "xdmf"
                )
                or source.xdmf.resolve().is_relative_to(directory)
            ):
                raise ValueError("Lobulus cache must be separate from its geometry")
        else:
            _, output_root = validate_roots(source.root, output_root)
        message, row = _build(source, output_root, profile, legacy, force)
        return ProcessingResult(
            source.study,
            source.run_id,
            "cached" if message == "cached" else "processed",
            message,
            time.perf_counter() - start,
            row,
        )
    except Exception as error:
        return ProcessingResult(
            source.study,
            source.run_id,
            "failed",
            f"{type(error).__name__}: {error}",
            time.perf_counter() - start,
            {},
        )


def _current_generation(source: RunSource, output_root: Path) -> Path:
    manifest = _manifest(source, output_root)
    if manifest is None:
        raise FileNotFoundError(f"No completed cache for {source.run_id}")
    generation = _generation(run_directory(source, output_root), manifest)
    if not _unchanged(manifest.get("sources"), source.root) or not _unchanged(
        manifest.get("outputs"), generation
    ):
        raise ValueError("Cache is stale; rerun processing")
    configuration = manifest.get("configuration")
    condition = asdict(source.condition) if source.condition is not None else None
    if (
        not isinstance(configuration, dict)
        or configuration.get("version") != CACHE_VERSION
        or configuration.get("condition") != condition
    ):
        raise ValueError("Cache schema changed; rerun processing")
    return generation


def load_processed_run(
    source: RunSource, output_root: Path
) -> tuple[xr.Dataset, xr.Dataset, xr.Dataset]:
    """Load completed cell, point, and summary data while holding the run lock."""
    if _manifest(source, output_root) is None:
        raise FileNotFoundError(f"No completed cache for {source.run_id}")
    with _run_lock(run_directory(source, output_root)):
        generation = _current_generation(source, output_root)
        return (
            xr.load_dataset(generation / "cells.nc", engine="scipy"),
            xr.load_dataset(generation / "points.nc", engine="scipy"),
            xr.load_dataset(generation / "summary.nc", engine="scipy"),
        )


def load_summary(source: RunSource, output_root: Path) -> xr.Dataset:
    """Load only compact native-time statistics without opening per-cell arrays."""
    if _manifest(source, output_root) is None:
        raise FileNotFoundError(f"No completed cache for {source.run_id}")
    with _run_lock(run_directory(source, output_root)):
        generation = _current_generation(source, output_root)
        return xr.load_dataset(generation / "summary.nc", engine="scipy")
