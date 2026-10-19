"""Serialize native summaries, matched-time statistics and spatial SPT profiles."""

from __future__ import annotations

import json
import tempfile
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import xarray as xr

from research.spt.artifacts import current_manifest, publish
from research.spt.cache import (
    SUMMARY_FIELDS,
    _current_generation,
    _generation,
    _run_lock,
    fingerprint,
    load_summary,
    open_processed_cells,
    run_directory,
)
from research.spt.plot_data import PlotSource, sample_times
from research.spt.sampling import blend, bracket

ANALYSIS_VERSION = 2
VOLUME = "element_volume_point_TPM"
POSITION = "rr_position"
DIAGNOSTICS = (VOLUME, "rr_Vext", "rr_Vli", "pressure")


@dataclass(frozen=True)
class AnalysisOptions:
    """Choose reproducible reductions independently of plot styling."""

    fields: tuple[str, ...] = SUMMARY_FIELDS
    time_points: int = 51
    bins: int = 20
    spatial_times: tuple[float, ...] = ()

    def validate(self) -> None:
        """Reject unsupported fields, duplicate times and invalid reduction sizes."""
        if (
            not self.fields
            or len(set(self.fields)) != len(self.fields)
            or set(self.fields) - set(SUMMARY_FIELDS)
        ):
            raise ValueError("Select unique SPT analysis scalar fields")
        if self.time_points < 2 or self.bins < 1:
            raise ValueError("time_points must be at least two; bins must be positive")
        if self.spatial_times and (
            not np.isfinite(self.spatial_times).all()
            or np.any(np.diff(self.spatial_times) <= 0)
        ):
            raise ValueError("Spatial times must be finite and strictly increasing")


def timelines(
    sources: list[PlotSource], options: AnalysisOptions
) -> dict[str, tuple[tuple[float, ...], tuple[float, ...]]]:
    """Keep study windows separate and retain each study's historical reference time."""
    options.validate()
    chosen = sample_times(sources, options.time_points)
    result = {}
    for family, (times, _) in chosen.items():
        start, end = times[0], times[-1]
        if start == end:
            raise ValueError("Analysis requires a nonzero common time interval")
        reference = 5400.0 if family == "convergence" else 36000.0
        middle = reference if start < reference < end else (start + end) / 2
        spatial = options.spatial_times or (start, middle, end)
        if spatial[0] < start or spatial[-1] > end:
            raise ValueError(f"Spatial times outside common {family} interval")
        result[family] = times, spatial
    return result


class CellSampler:
    """Keep two requested native cell frames from the local NetCDF cache."""

    def __init__(self, cells: xr.Dataset, fields: tuple[str, ...]) -> None:
        """Validate the cached coordinates and selected scalar dimensions."""
        self.cells = cells
        self.fields = fields
        self.times = np.asarray(cells.time.values, dtype=float)
        if (
            not len(self.times)
            or not np.isfinite(self.times).all()
            or np.any(np.diff(self.times) <= 0)
        ):
            raise ValueError("Invalid native time axis")
        missing = set(fields) - set(cells.data_vars)
        if missing:
            raise ValueError(
                f"Cache missing fields: {', '.join(sorted(missing))}; rerun cache processing"
            )
        if any(cells[name].dims != ("time", "cell") for name in fields):
            raise ValueError("Expected scalar time/cell arrays")
        self.frames: dict[int, dict[str, np.ndarray]] = {}
        self.frame_reads = 0

    def sample(self, time: float) -> dict[str, np.ndarray]:
        """Interpolate cells and volumes before any spatial reduction."""
        lower, upper, weight = bracket(self.times, time)
        self.frames = {
            key: value for key, value in self.frames.items() if key in {lower, upper}
        }
        for index in (lower, upper):
            if index not in self.frames:
                self.frames[index] = {
                    name: np.array(self.cells[name].isel(time=index).values, copy=True)
                    for name in self.fields
                }
                self.frame_reads += 1
        result = {
            name: blend(
                name, self.frames[lower][name], self.frames[upper][name], weight
            )
            for name in self.fields
        }
        if any(not np.isfinite(values).all() for values in result.values()):
            raise ValueError(f"Nonfinite selected cell values at {time:g} s")
        if "rr_necrosis" in result and (
            (result["rr_necrosis"] < 0).any() or (result["rr_necrosis"] > 1).any()
        ):
            raise ValueError("Necrosis values must lie between zero and one")
        return result


def weighted_statistics(
    values: np.ndarray, volumes: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Calculate volume-weighted population mean and standard deviation in float64."""
    volumes = np.asarray(volumes, dtype=np.float64)
    values = np.asarray(values, dtype=np.float64)
    if not np.isfinite(volumes).all() or (volumes < 0).any() or volumes.sum() <= 0:
        raise ValueError(
            "Cell volumes must be finite, nonnegative, with a positive total"
        )
    if not np.isfinite(values).all():
        raise ValueError("Analysis values must be finite")
    mean = np.sum(values * volumes, axis=-1) / volumes.sum()
    variance = (
        np.sum((values - mean[..., None]) ** 2 * volumes, axis=-1) / volumes.sum()
    )
    return mean, np.sqrt(variance)


def spatial_profile(
    frame: dict[str, np.ndarray], fields: tuple[str, ...], bins: int
) -> xr.Dataset:
    """Bin normalized source position, preserving empty bins as missing values."""
    position = frame[POSITION]
    if (
        not np.isfinite(position).all()
        or (position < -1e-6).any()
        or (position > 1 + 1e-6).any()
    ):
        raise ValueError("Normalized rr_position must lie between zero and one")
    volumes = frame[VOLUME]
    values = np.stack([frame[name] for name in fields])
    weighted_statistics(values, volumes)
    indices = np.minimum((np.clip(position, 0, 1) * bins).astype(int), bins - 1)
    means = np.full((len(fields), bins), np.nan)
    stds = means.copy()
    minima = means.copy()
    maxima = means.copy()
    count = np.bincount(indices, minlength=bins)
    totals = np.bincount(indices, weights=volumes.astype(float), minlength=bins)
    for index in np.flatnonzero(totals > 0):
        selected = indices == index
        means[:, index], stds[:, index] = weighted_statistics(
            values[:, selected], volumes[selected]
        )
        minima[:, index] = values[:, selected].min(axis=1)
        maxima[:, index] = values[:, selected].max(axis=1)
    edges = np.linspace(0, 1, bins + 1)
    return xr.Dataset(
        {
            "mean": (("field", "position"), means),
            "std": (("field", "position"), stds),
            "minimum": (("field", "position"), minima),
            "maximum": (("field", "position"), maxima),
            "cell_count": ("position", count.astype(np.int32)),
            "volume": ("position", totals),
            "bin_left": ("position", edges[:-1]),
            "bin_right": ("position", edges[1:]),
        },
        coords={"field": list(fields), "position": (edges[:-1] + edges[1:]) / 2},
        attrs={
            "weighting": "cell volume; population standard deviation",
            "empty_bins": "NaN statistics and zero count/volume",
            "position_definition": "source rr_position; 0 periportal, 1 perivenous",
        },
    )


def input_paths(source: PlotSource, data_root: Path) -> tuple[Path, ...]:
    """Depend on the validated numerical cache and its stable manifest pointer."""
    generation = _current_generation(source.run, data_root)
    return (
        run_directory(source.run, data_root) / "manifest.json",
        generation / "cells.nc",
        generation / "summary.nc",
    )


def analysis_configuration(
    source: PlotSource,
    options: AnalysisOptions,
    times: tuple[float, ...],
    spatial: tuple[float, ...],
) -> dict[str, object]:
    """Record mathematical choices and run conditions, without plot settings."""
    return {
        "version": ANALYSIS_VERSION,
        "source": str(source.run.xdmf),
        "condition": asdict(source.run.condition) if source.run.condition else None,
        "options": json.loads(json.dumps(asdict(options))),
        "times": list(times),
        "spatial_times": list(spatial),
    }


def prepare_analysis(
    source: PlotSource,
    data_root: Path,
    output_root: Path,
    options: AnalysisOptions,
    times: tuple[float, ...],
    spatial: tuple[float, ...],
    *,
    force: bool = False,
    dry_run: bool = False,
) -> dict[str, object] | None:
    """Publish compact reusable analysis data without reading network HDF arrays."""
    options.validate()
    directory = source.directory(output_root, "analysis-data")
    configuration = analysis_configuration(source, options, times, spatial)
    paths = input_paths(source, data_root)
    cached = current_manifest(directory, configuration, data_root)
    if dry_run or (cached and not force):
        return cached
    with _run_lock(directory):
        cached = current_manifest(directory, configuration, data_root)
        if cached and not force:
            return cached
        before = [fingerprint(path) for path in paths]
        all_fields = tuple(
            dict.fromkeys((*options.fields, POSITION, "rr_protein", *DIAGNOSTICS))
        )
        with load_summary(source.run, data_root) as native:
            required = [
                f"{field}_{stat}"
                for field in options.fields
                for stat in ("mean", "std")
            ]
            if set(required) - set(native.data_vars):
                raise ValueError(
                    "Missing native summary fields; rerun cache processing"
                )
            native = native[required].copy(deep=True)
        means, stds, volumes, profiles, frames = [], [], [], [], []
        selected_times = np.union1d(np.union1d(times, spatial), [source.start])
        initial = None
        with open_processed_cells(source.run, data_root) as cells:
            sampler = CellSampler(cells, all_fields)
            for time in selected_times:
                frame = sampler.sample(float(time))
                if time == source.start:
                    initial = spatial_profile(
                        frame, all_fields, options.bins
                    ).expand_dims(time=[source.start])
                if time in times:
                    mean, std = weighted_statistics(
                        np.stack([frame[name] for name in options.fields]),
                        frame[VOLUME],
                    )
                    means.append(mean)
                    stds.append(std)
                    volumes.append(np.sum(frame[VOLUME], dtype=float))
                if time in spatial:
                    profiles.append(spatial_profile(frame, all_fields, options.bins))
                    frames.append(
                        xr.Dataset(
                            {
                                name: ("cell", values.copy())
                                for name, values in frame.items()
                            },
                            coords={"cell": cells.cell.values.copy()},
                        )
                    )
            reads = sampler.frame_reads
        aligned = xr.Dataset(
            {
                "mean": (("time", "field"), np.stack(means)),
                "std": (("time", "field"), np.stack(stds)),
                "volume": ("time", volumes),
            },
            coords={"time": list(times), "field": list(options.fields)},
            attrs={
                "definition": "interpolate continuous cell values and volumes first, then reduce; hold previous discrete values; no extrapolation"
            },
        )
        profile = xr.concat(profiles, dim=xr.IndexVariable("time", list(spatial)))
        profile_cells = xr.concat(frames, dim=xr.IndexVariable("time", list(spatial)))
        assert initial is not None
        metadata: dict[str, object] = {
            "data_root": str(data_root),
            "title": source.title,
            "num_cells": len(profile_cells.cell),
            "native_frames_read": reads,
            "spatial_times": list(spatial),
            "source_start_s": source.start,
            "source_end_s": source.end,
        }
        with tempfile.TemporaryDirectory(
            prefix=".partial-", dir=directory
        ) as temporary:
            staging = Path(temporary)
            for name, dataset in (
                ("native", native),
                ("aligned", aligned),
                ("profiles", profile),
                ("cells", profile_cells),
                ("initial", initial),
            ):
                dataset.attrs.update(
                    study=source.run.study,
                    run_id=source.run.run_id,
                    time_units="s",
                    scalar_units="source concentration units; necrosis fraction; volumes retain source units",
                )
                dataset.to_netcdf(staging / f"{name}.nc", engine="scipy")
                dataset.close()
            publish(directory, staging, configuration, paths, before, metadata)
        return current_manifest(directory, configuration, data_root)


def load_analysis(source: PlotSource, output_root: Path, name: str) -> xr.Dataset:
    """Load a serialized run result while protecting its current generation."""
    if name not in {"native", "aligned", "profiles", "cells", "initial"}:
        raise ValueError(f"Unknown analysis dataset: {name}")
    directory = source.directory(output_root, "analysis-data")
    if not (directory / "manifest.json").is_file():
        raise FileNotFoundError(f"No completed analysis for {source.run.run_id}")
    with _run_lock(directory):
        manifest = json.loads((directory / "manifest.json").read_text())
        configuration = manifest["configuration"]
        _current_generation(source.run, Path(manifest["metadata"]["data_root"]))
        if (
            configuration.get("version") != ANALYSIS_VERSION
            or current_manifest(
                directory, configuration, Path(manifest["metadata"]["data_root"])
            )
            is None
        ):
            raise ValueError("Analysis is stale; rerun analysis")
        generation = _generation(directory, manifest)
        return xr.load_dataset(generation / f"{name}.nc", engine="scipy")
