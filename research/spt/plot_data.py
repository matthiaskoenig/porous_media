"""Select SPT geometry sources and sample only requested scalar fields and times."""

from __future__ import annotations

import hashlib
import json
import tempfile
from contextlib import ExitStack
from dataclasses import dataclass
from pathlib import Path
from typing import Literal, cast

import h5py
import meshio
import numpy as np

from research.spt.artifacts import (
    atomic_json,
    current_manifest,
    fingerprint,
    publish,
    run_lock,
)
from research.spt.cache import _read_frame, cached_catalog, validate_roots
from research.spt.catalog import (
    CatalogValue,
    ParsedRun,
    RunSource,
    Study,
    discover_runs,
    inspect_run,
)
from research.spt.cli import selected_studies
from research.spt.conditions import pattern_idx2name
from research.spt.geometry import geometry_jobs, geometry_path
from research.spt.sampling import blend, bracket

SAMPLE_VERSION = 1
type Geometry = Literal["sixth", "lobulus", "combined"]
GEOMETRIES: tuple[Geometry, ...] = ("sixth", "lobulus", "combined")


@dataclass(frozen=True)
class PlotSource:
    """Describe a selected history without reading its numerical arrays."""

    run: RunSource
    geometry: Geometry
    start: float
    end: float
    title: str
    sector_labels: tuple[str, ...] = ()
    extra_sources: tuple[Path, ...] = ()
    dose: int | None = None

    @property
    def family(self) -> str:
        """Return the study family used for shared times and color scales."""
        return self.run.study.split("_")[0]

    def directory(self, output_root: Path, category: str) -> Path:
        """Keep sample caches and figures in separate safe local directories."""
        root = output_root.expanduser().resolve()
        path = root / category / self.run.study / self.run.run_id
        if not path.is_relative_to(root) or path.resolve() != path.absolute():
            raise ValueError("Plot directories must not traverse symbolic links")
        return path


def _time(row: dict[str, CatalogValue], name: str) -> float:
    value = row.get(name)
    if not isinstance(value, (int, float)) or not np.isfinite(value):
        raise ValueError(f"Invalid catalogue time: {name}")
    return float(value)


def _title(run: RunSource, geometry: Geometry) -> str:
    condition = run.condition
    if condition is None:
        return f"{run.run_id} | {geometry}"
    return f"{run.run_id} | {pattern_idx2name[condition.pattern_key]} | dose {condition.boundary_flow_key + 1} | {geometry}"


def source_row(
    source: RunSource, data_root: Path, plot_root: Path
) -> dict[str, CatalogValue]:
    """Reuse validated time bounds without reopening unchanged descriptors."""
    row = cached_catalog(source, data_root)
    if row:
        return row
    directory = plot_root / "plot-data" / source.study / source.run_id
    try:
        stored = json.loads((directory / "manifest.json").read_text())
        configuration = stored["configuration"]
        if configuration["version"] == SAMPLE_VERSION and configuration[
            "descriptor"
        ] == str(source.xdmf):
            current = current_manifest(directory, configuration, source.root)
            if current:
                metadata = cast(dict[str, object], current["metadata"])
                start = metadata.get("source_start_s")
                end = metadata.get("source_end_s")
                if isinstance(start, (int, float)) and isinstance(end, (int, float)):
                    return {"tstart_s": float(start), "tend_s": float(end)}
    except (OSError, ValueError, KeyError, TypeError):
        pass
    return inspect_run(source).row


def plot_sources(
    input_root: Path,
    data_root: Path,
    *,
    study: str = "both",
    geometry: str = "all",
    run_ids: tuple[str, ...] = (),
    plot_root: Path | None = None,
) -> list[PlotSource]:
    """Resolve raw sixths and published local geometries using stable condition IDs."""
    input_root, data_root = validate_roots(input_root, data_root)
    if geometry not in {"all", *GEOMETRIES}:
        raise ValueError(f"Unknown geometry: {geometry}")
    kinds = GEOMETRIES if geometry == "all" else (geometry,)
    raw = discover_runs(input_root, selected_studies(study))
    selected_raw = [source for source in raw if not run_ids or source.run_id in run_ids]
    result = []
    if "sixth" in kinds:
        for source in selected_raw:
            row = source_row(source, data_root, plot_root or data_root)
            result.append(
                PlotSource(
                    source,
                    "sixth",
                    _time(row, "tstart_s"),
                    _time(row, "tend_s"),
                    _title(source, "sixth"),
                    dose=source.condition.boundary_flow_key
                    if source.condition
                    else None,
                )
            )
    if "lobulus" in kinds:
        for source in selected_raw:
            folder = (
                data_root
                / source.study.replace("_sixth", "_lobulus")
                / "xdmf"
                / source.run_id
            )
            path = geometry_path(folder)
            run = RunSource(
                data_root,
                cast(Study, source.study.replace("_sixth", "_lobulus")),
                source.run_id,
                source.condition,
                source.condition_source,
                path,
                path.parent / "model.feb",
            )
            manifest = json.loads((folder / "manifest.json").read_text())
            row = manifest["catalog"]
            result.append(
                PlotSource(
                    run,
                    "lobulus",
                    _time(row, "tstart_s"),
                    _time(row, "tend_s"),
                    _title(run, "lobulus"),
                    extra_sources=(folder / "manifest.json",),
                    dose=source.condition.boundary_flow_key
                    if source.condition
                    else None,
                )
            )
    if "combined" in kinds and (
        not run_ids
        or any(
            name.startswith("dose_") or name == "convergence_mixed" for name in run_ids
        )
    ):
        for job in geometry_jobs(input_root, study=study, lobulus=False):
            if run_ids and job.run_id not in run_ids:
                continue
            folder = job.directory(data_root)
            path = geometry_path(folder)
            manifest = json.loads((folder / "manifest.json").read_text())
            row = manifest["catalog"]
            condition = job.sources[0].condition
            labels = tuple(
                pattern_idx2name[source.condition.pattern_key]
                if source.condition
                else source.run_id.removeprefix("convergence_")
                for source in job.sources
            )
            title = (
                f"{job.run_id} | six zonation patterns"
                if condition
                else "Convergence | six resolutions"
            )
            run = RunSource(
                data_root,
                cast(Study, job.study),
                job.run_id,
                None,
                "published geometry",
                path,
                path.parent / "model.feb",
            )
            result.append(
                PlotSource(
                    run,
                    "combined",
                    _time(row, "tstart_s"),
                    _time(row, "tend_s"),
                    title,
                    labels,
                    (folder / "manifest.json",),
                    condition.boundary_flow_key if condition else None,
                )
            )
    unknown = set(run_ids) - {source.run.run_id for source in result}
    if unknown or not result:
        raise ValueError(
            f"No selected plot histories or unknown runs: {', '.join(sorted(unknown))}"
        )
    return sorted(result, key=lambda source: (source.run.study, source.run.run_id))


def sample_times(
    sources: list[PlotSource],
    snapshots: int,
    explicit: tuple[float, ...] = (),
    video_frames: int = 0,
) -> dict[str, tuple[tuple[float, ...], tuple[float, ...]]]:
    """Choose common ranges separately for simulation and convergence histories."""
    if snapshots < 1 or video_frames < 0 or video_frames == 1:
        raise ValueError(
            "snapshots must be positive; video_frames must be zero or at least two"
        )
    if explicit and (not np.isfinite(explicit).all() or np.any(np.diff(explicit) <= 0)):
        raise ValueError("Snapshot times must be finite and strictly increasing")
    result = {}
    for family in {source.family for source in sources}:
        selected = [source for source in sources if source.family == family]
        start, end = (
            max(source.start for source in selected),
            min(source.end for source in selected),
        )
        if start > end:
            raise ValueError(f"No common time interval for {family}")
        times = explicit or (
            (end,)
            if snapshots == 1
            else tuple(float(value) for value in np.linspace(start, end, snapshots))
        )
        if min(times) < start or max(times) > end:
            raise ValueError(
                f"Requested times outside the common {family} range [{start}, {end}] s"
            )
        video = (
            tuple(float(value) for value in np.linspace(start, end, video_frames))
            if video_frames
            else ()
        )
        result[family] = (tuple(dict.fromkeys(times)), video)
    return result


class ScalarSampler:
    """Read only selected cell scalars with at most two native frames in memory."""

    def __init__(
        self, parsed: ParsedRun, fields: tuple[str, ...], handles: dict[str, h5py.File]
    ) -> None:
        """Validate the selected scalar schema without reading unused fields."""
        self.parsed = parsed
        self.fields = fields
        self.handles = handles
        self.frames: dict[int, dict[str, np.ndarray]] = {}
        self.field_reads = 0
        missing = set(fields) - set(parsed.fields["Cell"])
        if missing:
            raise ValueError(
                f"Missing scalar cell fields: {', '.join(sorted(missing))}"
            )

    def sample(self, time: float) -> np.ndarray:
        """Interpolate continuous fields and hold discrete values without extrapolation."""
        lower, upper, weight = bracket(self.parsed.times, time)
        self.frames = {
            index: data
            for index, data in self.frames.items()
            if index in {lower, upper}
        }
        for index in {lower, upper}:
            if index not in self.frames:
                self.frames[index] = _read_frame(
                    self.parsed, index, "Cell", self.fields, self.handles
                )
                self.field_reads += len(self.fields)
        arrays = []
        for field in self.fields:
            low = self.frames[lower][field]
            values = blend(field, low, self.frames[upper][field], weight)
            if (
                values.shape != (self.parsed.count("Cell"),)
                or not np.isfinite(values).all()
            ):
                raise ValueError(f"Invalid scalar plot values: {field} at {time:g} s")
            arrays.append(values)
        return np.stack(arrays)


def projected_faces(
    points: np.ndarray, cells: list[meshio.CellBlock]
) -> tuple[np.ndarray, np.ndarray]:
    """Triangulate planar faces while preserving their original cell identities."""
    if points.ndim != 2 or points.shape[1] != 3 or not np.isfinite(points).all():
        raise ValueError("Plot geometry requires finite three-dimensional coordinates")
    tolerance = max(float(np.ptp(points, axis=0).max()) * 1e-7, np.finfo(float).eps)
    faces = []
    ids = []
    offset = 0
    for block in cells:
        vertices = points[block.data]
        if block.type in {"wedge", "hexahedron"}:
            count = 3 if block.type == "wedge" else 4
            if not np.allclose(
                vertices[:, :count, :2], vertices[:, count:, :2], atol=tolerance, rtol=0
            ):
                raise ValueError("Planar plotting requires vertical mesh extrusions")
            base = vertices[:, :count]
        elif block.type in {"triangle", "quad"}:
            count = 3 if block.type == "triangle" else 4
            base = vertices
        else:
            raise ValueError(
                f"Unsupported planar cell type: {block.type}; use the library's 3D renderer"
            )
        if np.ptp(base[:, :, 2]) > tolerance:
            raise ValueError("Plot faces must lie in one XY plane")
        for indices in ((0, 1, 2),) if count == 3 else ((0, 1, 2), (0, 2, 3)):
            faces.append(base[:, indices, :2])
            ids.append(np.arange(offset, offset + len(block), dtype=np.int64))
        offset += len(block)
    if not faces:
        raise ValueError("No plot faces found")
    return np.concatenate(faces), np.concatenate(ids)


def prepare_samples(
    source: PlotSource,
    output_root: Path,
    fields: tuple[str, ...],
    times: tuple[float, ...],
    *,
    force: bool = False,
    dry_run: bool = False,
) -> dict[str, object] | None:
    """Cache compact requested samples; never write interpolated simulation histories."""
    if not fields or len(set(fields)) != len(fields) or not times:
        raise ValueError("Fields and times must be nonempty and unique")
    directory = source.directory(output_root, "plot-data")
    configuration: dict[str, object] = {
        "version": SAMPLE_VERSION,
        "descriptor": str(source.run.xdmf),
        "fields": list(fields),
        "times": list(times),
    }
    cached = current_manifest(directory, configuration, source.run.root)
    if dry_run or (cached and not force):
        return cached
    with run_lock(directory):
        cached = current_manifest(directory, configuration, source.run.root)
        if cached and not force:
            return cached
        parsed = inspect_run(source.run)
        paths = tuple(sorted({*parsed.sources, *source.extra_sources}))
        before = [fingerprint(path) for path in paths]
        with meshio.xdmf.TimeSeriesReader(source.run.xdmf) as reader:
            points, cells = reader.read_points_cells()
        faces, cell_ids = projected_faces(points, cells)
        with ExitStack() as stack:
            handles = {
                name: stack.enter_context(h5py.File(path, "r"))
                for name, path in parsed.references.items()
            }
            sampler = ScalarSampler(parsed, fields, handles)
            values = np.stack([sampler.sample(time) for time in times])
            sectors = (
                _read_frame(parsed, 0, "Cell", ("sector_id",), handles)["sector_id"]
                if "sector_id" in parsed.fields["Cell"]
                else np.full(parsed.count("Cell"), -1, dtype=np.int32)
            )
        metadata: dict[str, object] = {
            "title": source.title,
            "fields": list(fields),
            "times": list(times),
            "source_start_s": float(parsed.times[0]),
            "source_end_s": float(parsed.times[-1]),
            "num_cells": parsed.count("Cell"),
            "num_faces": len(faces),
            "mesh_hash": hashlib.blake2b(
                faces.tobytes() + cell_ids.tobytes(), digest_size=16
            ).hexdigest(),
            "ranges": {
                name: [float(values[:, index].min()), float(values[:, index].max())]
                for index, name in enumerate(fields)
            },
            "sampling": "linear continuous; previous discrete; no extrapolation",
            "hdf_field_reads": sampler.field_reads,
        }
        if before != [fingerprint(path) for path in paths]:
            raise RuntimeError("Inputs changed while sampling plot data")
        with tempfile.TemporaryDirectory(
            prefix=".partial-", dir=directory
        ) as temporary:
            staging = Path(temporary)
            np.savez(
                staging / "samples.npz",
                faces=faces,
                cell_ids=cell_ids,
                values=values,
                sectors=sectors,
                times=np.array(times),
            )
            atomic_json(staging / "metadata.json", metadata)
            if before != [fingerprint(path) for path in paths]:
                raise RuntimeError("Inputs changed before plot sample publication")
            publish(directory, staging, configuration, paths, before, metadata)
        return current_manifest(directory, configuration, source.run.root)
