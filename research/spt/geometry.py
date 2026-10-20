"""Reconstruct lobules and independent sector comparisons from SPT time series.

Derived XDMF/HDF pairs are compressed, self-contained and published through an
atomic manifest. Mixed comparisons are visualization datasets, not new solver
runs or inputs to aggregate quantitative analysis.
"""

from __future__ import annotations

import hashlib
import json
import shutil
import tempfile
import time
from contextlib import ExitStack
from dataclasses import asdict, dataclass
from functools import partial
from pathlib import Path
from typing import cast

import meshio
import numpy as np
import pandas as pd
from numpy.typing import ArrayLike

from porous_media.data.xdmf_tools import XDMFTimeSeriesWriter
from research.geometry_conversion import (
    find_reflection_angle_in_sixth,
    sector_transforms,
    transform_cells,
    transform_fields,
)
from research.spt.artifacts import (
    atomic_json,
    fingerprint,
    generation_path,
    publish,
    records_current,
    run_lock,
)
from research.spt.cache import ProcessingResult, validate_roots
from research.spt.catalog import (
    CatalogValue,
    RunSource,
    Study,
    discover_runs,
    inspect_run,
)
from research.spt.cli import selected_studies, write_catalog
from research.spt.conditions import conditions_from_source
from research.spt.execution import execute_jobs, print_result
from research.spt.sampling import DISCRETE_FIELDS as DISCRETE_FIELDS
from research.spt.sampling import blend, bracket

GEOMETRY_VERSION = 1
CONVERGENCE_MESHES = ("06250", "09375", "12500", "18750", "25000", "37250")


@dataclass(frozen=True)
class GeometryJob:
    """Identify six ordered source sectors and the purpose of their composition."""

    study: str
    run_id: str
    sources: tuple[RunSource, ...]
    mixed: bool = False

    def directory(self, output_root: Path) -> Path:
        """Return a stable manifest directory without following output links."""
        output_root = output_root.expanduser().resolve()
        directory = output_root / self.study / "xdmf" / self.run_id
        if (
            not directory.is_relative_to(output_root)
            or directory.resolve() != directory.absolute()
        ):
            raise ValueError("Geometry directories must not traverse symbolic links")
        return directory

    def configuration(self) -> dict[str, object]:
        """Record choices that affect geometry, fields or temporal alignment."""
        return {
            "version": GEOMETRY_VERSION,
            "sources": [str(source.xdmf) for source in self.sources],
            "conditions": [
                asdict(source.condition) if source.condition else None
                for source in self.sources
            ],
            "mixed": self.mixed,
            "time_policy": "first sector native times within common interval; linear continuous, previous discrete"
            if self.mixed
            else "native",
            "seams": "separate points; alternating reflections; source orientation preserved",
            "compression": "gzip level 1 with shuffle; identical arrays share HDF storage",
        }


class _CompressedWriter(XDMFTimeSeriesWriter):
    """Reuse identical arrays and compress HDF datasets within one local sidecar."""

    def __enter__(self) -> _CompressedWriter:
        super().__enter__()
        self._datasets: dict[tuple[str, tuple[int, ...], bytes], str] = {}
        return self

    def numpy_to_xml_string(self, data: np.ndarray) -> str:
        """Write meshio arrays with portable HDF compression and hard links."""
        data = np.ascontiguousarray(data)
        name = f"data{self.data_counter}"
        self.data_counter += 1
        key = (
            data.dtype.str,
            data.shape,
            hashlib.blake2b(data, digest_size=32).digest(),
        )
        previous = self._datasets.get(key)
        if previous is None:
            self.h5_file.create_dataset(
                name, data=data, compression="gzip", compression_opts=1, shuffle=True
            )
            self._datasets[key] = name
        else:
            self.h5_file[name] = self.h5_file[previous]
        return f"{Path(self.h5_filename).name}:/{name}"


def geometry_jobs(
    input_root: Path,
    *,
    study: str = "both",
    run_ids: tuple[str, ...] = (),
    lobulus: bool = True,
    combinations: bool = True,
) -> list[GeometryJob]:
    """Select full lobules and six-pattern views using archived conditions."""
    runs = discover_runs(input_root, selected_studies(study), run_ids)
    if not run_ids and study in {"both", "simulations"}:
        conditions, _ = conditions_from_source(input_root)
        missing = {condition.sim_id for condition in conditions} - {
            source.run_id for source in runs
        }
        if missing:
            raise ValueError(
                f"Missing expected simulations: {', '.join(sorted(missing))}"
            )
    jobs = []
    if lobulus:
        jobs.extend(
            GeometryJob(
                source.study.replace("_sixth", "_lobulus"), source.run_id, (source,) * 6
            )
            for source in runs
        )
    if not combinations:
        return jobs
    if run_ids:
        raise ValueError("Combinations require complete studies; omit --runs")
    by_id = {source.run_id: source for source in runs}
    if study in {"both", "convergence"}:
        names = [f"convergence_{key}" for key in CONVERGENCE_MESHES]
        missing = set(names) - by_id.keys()
        if missing:
            raise ValueError(
                f"Missing convergence comparison meshes: {', '.join(sorted(missing))}"
            )
        jobs.append(
            GeometryJob(
                "convergence_combinations",
                "convergence_mixed",
                tuple(by_id[name] for name in names),
                True,
            )
        )
    if study in {"both", "simulations"}:
        conditions, _ = conditions_from_source(input_root)
        missing = {condition.sim_id for condition in conditions} - by_id.keys()
        if missing:
            raise ValueError(
                f"Missing expected simulations: {', '.join(sorted(missing))}"
            )
        for dose in sorted({condition.boundary_flow_key for condition in conditions}):
            selected = sorted(
                (
                    source
                    for source in runs
                    if source.condition and source.condition.boundary_flow_key == dose
                ),
                key=lambda source: (
                    source.condition.pattern_key if source.condition else -1
                ),
            )
            if len(selected) != 6 or {
                source.condition.pattern_key for source in selected if source.condition
            } != set(range(6)):
                raise ValueError("Each dose requires all six zonation patterns")
            jobs.append(
                GeometryJob(
                    "simulations_combinations",
                    f"dose_{dose + 1:02}",
                    tuple(selected),
                    True,
                )
            )
    return jobs


def geometry_manifest(directory: Path) -> dict[str, object]:
    """Require a complete, current geometry generation and unchanged sources."""
    directory = directory.expanduser().absolute()
    if directory.resolve() != directory:
        raise ValueError("Geometry directories must not traverse symbolic links")
    manifest = json.loads((directory / "manifest.json").read_text(encoding="utf-8"))
    if not isinstance(manifest, dict):
        raise ValueError("Invalid geometry manifest")
    configuration = manifest.get("configuration")
    if (
        not isinstance(configuration, dict)
        or configuration.get("version") != GEOMETRY_VERSION
    ):
        raise ValueError("Geometry schema changed; rerun geometry processing")
    path = generation_path(directory, manifest)
    if not records_current(manifest.get("sources")) or not records_current(
        manifest.get("outputs"), root=path
    ):
        raise ValueError("Geometry is stale; rerun geometry processing")
    return manifest


def geometry_path(directory: Path) -> Path:
    """Resolve a published descriptor for analysis or visualization."""
    directory = directory.expanduser().absolute()
    manifest = geometry_manifest(directory)
    generation = directory / str(manifest["generation"])
    paths = list(generation.glob("*.xdmf"))
    if len(paths) != 1:
        raise ValueError("Expected one published XDMF descriptor")
    return paths[0]


def _status(job: GeometryJob, output_root: Path) -> dict[str, object] | None:
    try:
        manifest = geometry_manifest(job.directory(output_root))
    except (OSError, ValueError, KeyError):
        return None
    return manifest if manifest.get("configuration") == job.configuration() else None


class _Sampler:
    """Keep at most two source frames for interpolation on an explicit time axis."""

    def __init__(self, reader: meshio.xdmf.TimeSeriesReader, times: np.ndarray) -> None:
        self.reader = reader
        self.times = times
        self.frames: dict[
            int, tuple[dict[str, np.ndarray], dict[str, list[np.ndarray]]]
        ] = {}

    def sample(
        self, time_value: float
    ) -> tuple[dict[str, np.ndarray], dict[str, list[np.ndarray]]]:
        """Interpolate continuous fields and hold discrete fields without extrapolation."""
        lower, upper, weight = bracket(self.times, time_value)
        self.frames = {
            key: value for key, value in self.frames.items() if key in {lower, upper}
        }
        for step in (lower, upper):
            if step not in self.frames:
                _, points, cells = self.reader.read_data(step)
                self.frames[step] = (points, cells)
        points_low, cells_low = self.frames[lower]
        if upper == lower:
            return points_low, cells_low
        points_up, cells_up = self.frames[upper]
        if points_low.keys() != points_up.keys() or cells_low.keys() != cells_up.keys():
            raise ValueError("Field schema changes between time steps")
        return (
            {
                name: blend(name, values, points_up[name], weight)
                for name, values in points_low.items()
            },
            {
                name: [
                    blend(name, low, up, weight)
                    for low, up in zip(blocks, cells_up[name], strict=True)
                ]
                for name, blocks in cells_low.items()
            },
        )


def write_composition(
    sources: tuple[RunSource, ...], output: Path
) -> dict[str, object]:
    """Write six independent sectors, preserving all source fields and their types.

    Mixed runs use the first sector's native times within the common interval,
    including its endpoints when needed. Continuous fields interpolate linearly;
    categorical fields and necrosis flags hold their previous sampled state.
    """
    if len(sources) != 6:
        raise ValueError("A lobulus requires six source sectors")
    output = output.expanduser().absolute()
    if any(
        output.resolve().is_relative_to(source.root.resolve()) for source in sources
    ):
        raise ValueError("Geometry output must be outside the input tree")
    if any(
        path.exists() or path.is_symlink()
        for path in (output, output.with_suffix(".h5"))
    ):
        raise FileExistsError(f"Geometry output already exists: {output}")
    unique = tuple(dict.fromkeys(sources))
    parsed = {source: inspect_run(source) for source in unique}
    start = max(run.times[0] for run in parsed.values())
    end = min(run.times[-1] for run in parsed.values())
    if start > end:
        raise ValueError("Source time ranges do not overlap")
    reference = parsed[sources[0]].times
    times = np.unique(
        np.concatenate(
            ([start], reference[(reference >= start) & (reference <= end)], [end])
        )
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    limits: dict[str, tuple[float, float]] = {}
    with ExitStack() as stack:
        # Sequential copies avoid thousands of small HDF reads over a network
        # filesystem. Keep the relative source layout so XML needs no rewriting.
        local = Path(
            stack.enter_context(
                tempfile.TemporaryDirectory(prefix=".sources-", dir=output.parent)
            )
        )
        roots = {
            root: local / f"input-{index}"
            for index, root in enumerate(
                dict.fromkeys(source.root for source in unique)
            )
        }
        copied: set[Path] = set()
        for source, run in parsed.items():
            for path in (source.xdmf, *run.references.values()):
                target = roots[source.root] / path.relative_to(source.root)
                if target not in copied:
                    target.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copyfile(path, target)
                    copied.add(target)
        readers = {
            source: stack.enter_context(
                meshio.xdmf.TimeSeriesReader(
                    roots[source.root] / source.xdmf.relative_to(source.root)
                )
            )
            for source in unique
        }
        meshes = {
            source: reader.read_points_cells() for source, reader in readers.items()
        }
        samplers = {
            source: _Sampler(readers[source], run.times)
            for source, run in parsed.items()
        }
        matrices = []
        points = []
        cells: list[meshio.CellBlock] = []
        sector_cells: list[ArrayLike] = []
        source_cells: list[ArrayLike] = []
        sector_points = []
        offset = 0
        for index, source in enumerate(sources):
            coordinates, blocks = meshes[source]
            angle = find_reflection_angle_in_sixth(coordinates)
            # Each archived mesh occupies the same 60-degree sector.
            radii = np.linalg.norm(coordinates[:, :2], axis=1)
            angles = np.angle(
                np.exp(1j * (np.arctan2(coordinates[:, 1], coordinates[:, 0]) - angle))
            )
            angles = angles[radii > radii.max() * 1e-8]
            if not np.isclose(np.ptp(angles), np.pi / 3, atol=2e-5):
                raise ValueError(
                    "SPT reconstruction requires a 60-degree source sector"
                )
            if matrices and not np.isclose(
                angle, find_reflection_angle_in_sixth(meshes[sources[0]][0]), atol=2e-5
            ):
                raise ValueError("Source sectors have different orientations")
            matrix = sector_transforms(angle)[index]
            matrices.append(matrix)
            points.append(coordinates @ matrix.T)
            sector_points.append(np.full(len(coordinates), index, dtype=np.int32))
            source_offset = 0
            for block in blocks:
                cells.append(
                    meshio.CellBlock(
                        block.type, transform_cells(block, matrix) + offset
                    )
                )
                sector_cells.append(np.full(len(block), index, dtype=np.int32))
                source_cells.append(
                    np.arange(source_offset, source_offset + len(block), dtype=np.int32)
                )
                source_offset += len(block)
            offset += len(coordinates)
        with _CompressedWriter(output, data_format="HDF") as writer:
            writer.write_points_cells(
                np.concatenate(points), [(block.type, block.data) for block in cells]
            )
            for time_value in times:
                frames = {
                    source: sampler.sample(float(time_value))
                    for source, sampler in samplers.items()
                }
                point_sectors = []
                cell_sectors: list[dict[str, np.ndarray]] = []
                for source, matrix in zip(sources, matrices, strict=True):
                    point_data, cell_data = frames[source]
                    if any(
                        name in point_data or name in cell_data
                        for name in ("sector_id", "source_cell_id")
                    ):
                        raise ValueError(
                            "Source fields conflict with composition metadata"
                        )
                    point_sectors.append(transform_fields(point_data, matrix))
                    for block_index in range(len(meshes[source][1])):
                        cell_sectors.append(
                            transform_fields(
                                {
                                    name: blocks[block_index]
                                    for name, blocks in cell_data.items()
                                },
                                matrix,
                            )
                        )
                if any(
                    sector.keys() != point_sectors[0].keys() for sector in point_sectors
                ) or any(
                    sector.keys() != cell_sectors[0].keys() for sector in cell_sectors
                ):
                    raise ValueError("Source sectors have different field schemas")
                composed_points: dict[str, ArrayLike] = {
                    name: np.concatenate([sector[name] for sector in point_sectors])
                    for name in point_sectors[0]
                }
                composed_cells: dict[str, list[ArrayLike]] = {
                    name: [sector[name] for sector in cell_sectors]
                    for name in cell_sectors[0]
                }
                composed_points["sector_id"] = np.concatenate(sector_points)
                composed_cells["sector_id"] = sector_cells
                composed_cells["source_cell_id"] = source_cells
                for fields in (
                    composed_points,
                    {
                        name: np.concatenate(blocks)
                        for name, blocks in composed_cells.items()
                    },
                ):
                    for name, values in fields.items():
                        array = np.asarray(values)
                        if not np.isfinite(array).all():
                            raise ValueError(f"Nonfinite geometry field: {name}")
                        previous = limits.get(name, (float("inf"), float("-inf")))
                        limits[name] = (
                            min(previous[0], float(array.min())),
                            max(previous[1], float(array.max())),
                        )
                writer.write_data(
                    float(time_value),
                    point_data=composed_points,
                    cell_data=composed_cells,
                )
    atomic_json(output.with_name(f"{output.stem}_limits.json"), {"limits": limits})
    return {
        "num_points": offset,
        "num_cells": sum(len(block) for block in cells),
        "num_steps": len(times),
        "tstart_s": float(start),
        "tend_s": float(end),
        "sector_sources": [source.run_id for source in sources],
        "source_end_times_s": [float(parsed[source].times[-1]) for source in sources],
    }


def _geometry_catalog(
    job: GeometryJob, metadata: dict[str, object], generation: Path
) -> dict[str, CatalogValue]:
    return {
        "study": job.study,
        "run_id": job.run_id,
        "xdmf": str(generation / f"{job.run_id}.xdmf"),
        "num_points": cast(int, metadata["num_points"]),
        "num_cells": cast(int, metadata["num_cells"]),
        "num_steps": cast(int, metadata["num_steps"]),
        "tstart_s": cast(float, metadata["tstart_s"]),
        "tend_s": cast(float, metadata["tend_s"]),
        "quantitative_analysis": not job.mixed,
    }


def process_geometry(
    job: GeometryJob, output_root: Path, force: bool = False
) -> ProcessingResult:
    """Build one atomic generation; a failed rebuild preserves its predecessor."""
    started = time.perf_counter()
    try:
        for source in job.sources:
            _, output_root = validate_roots(source.root, output_root)
        directory = job.directory(output_root)
        with run_lock(directory):
            current = _status(job, output_root)
            if current is not None and not force:
                return ProcessingResult(
                    job.study,
                    job.run_id,
                    "cached",
                    "cached",
                    time.perf_counter() - started,
                    cast(dict[str, CatalogValue], current["catalog"]),
                )
            files = {
                path
                for source in set(job.sources)
                for path in inspect_run(source).sources
            }
            before = [fingerprint(path) for path in sorted(files)]
            with tempfile.TemporaryDirectory(
                prefix=".partial-", dir=directory
            ) as temporary:
                staging = Path(temporary)
                metadata = write_composition(
                    job.sources, staging / f"{job.run_id}.xdmf"
                )
                preview_metadata: dict[str, object] | None = None
                if job.mixed:
                    from research.spt.geometry_figures import (
                        COMPARISON_PLOT_VERSION,
                        render_comparison,
                    )

                    render_comparison(job, staging / f"{job.run_id}.xdmf")
                    preview_metadata = {
                        "comparison_plot_version": COMPARISON_PLOT_VERSION
                    }
                material = job.sources[0].material
                if not job.mixed and material.is_file():
                    shutil.copyfile(material, staging / "model.feb")
                atomic_json(
                    staging / "provenance.json",
                    {
                        **job.configuration(),
                        **metadata,
                        "quantitative_analysis": not job.mixed,
                    },
                )
                if before != [fingerprint(path) for path in sorted(files)]:
                    raise RuntimeError("Inputs changed during geometry processing")
                generation = publish(
                    directory,
                    staging,
                    job.configuration(),
                    tuple(sorted(files)),
                    before,
                    preview_metadata,
                    catalog=lambda path: _geometry_catalog(job, metadata, path),
                )
                row = _geometry_catalog(job, metadata, generation)
        return ProcessingResult(
            job.study,
            job.run_id,
            "processed",
            "reconstructed",
            time.perf_counter() - started,
            row,
        )
    except Exception as error:
        return ProcessingResult(
            job.study,
            job.run_id,
            "failed",
            f"{type(error).__name__}: {error}",
            time.perf_counter() - started,
            {},
        )


def process_geometries(
    input_root: Path,
    output_root: Path,
    *,
    study: str = "both",
    run_ids: tuple[str, ...] = (),
    lobulus: bool = True,
    combinations: bool = True,
    workers: int = 2,
    force: bool = False,
    dry_run: bool = False,
) -> list[ProcessingResult]:
    """Build independent geometries with bounded workers and cheap resume checks."""
    if workers < 1:
        raise ValueError("workers must be positive")
    input_root, output_root = validate_roots(input_root, output_root)
    jobs = geometry_jobs(
        input_root,
        study=study,
        run_ids=run_ids,
        lobulus=lobulus,
        combinations=combinations,
    )
    results = []
    pending = []
    for job in jobs:
        manifest = _status(job, output_root)
        if manifest is not None and job.mixed and not force and not dry_run:
            from research.spt.geometry_figures import refresh_comparison

            try:
                manifest = refresh_comparison(job, output_root)
            except Exception as error:
                results.append(ProcessingResult.failed(job.study, job.run_id, error))
                continue
        if dry_run or (manifest is not None and not force):
            results.append(
                ProcessingResult(
                    job.study,
                    job.run_id,
                    "cached" if manifest is not None and not force else "pending",
                    "cached"
                    if manifest is not None and not force
                    else "reconstruction required",
                    0.0,
                    cast(dict[str, CatalogValue], manifest["catalog"])
                    if manifest
                    else {},
                )
            )
        else:
            pending.append(job)
    for _, result in execute_jobs(
        pending,
        partial(process_geometry, output_root=output_root, force=force),
        workers=workers,
        failed=lambda job, error: ProcessingResult.failed(job.study, job.run_id, error),
    ):
        results.append(result)
        print_result(result)
    results.sort(key=lambda result: (result.study, result.run_id))
    for result in results:
        if result.row and result.study.endswith("_combinations"):
            descriptor = Path(str(result.row["xdmf"]))
            result.row["figure_png"] = str(descriptor.with_suffix(".png"))
            result.row["figure_pdf"] = str(descriptor.with_suffix(".pdf"))
    if not dry_run:
        atomic_json(
            output_root / "geometry-report.json", [asdict(result) for result in results]
        )
        write_catalog(
            pd.DataFrame([result.row for result in results if result.row]),
            output_root,
            name="geometry-catalog",
        )
    return results


def lobulus_sources(
    input_root: Path,
    output_root: Path,
    *,
    study: str = "both",
    run_ids: tuple[str, ...] = (),
) -> list[RunSource]:
    """Resolve homogeneous lobules for quantitative analysis, excluding mixed views."""
    _, output_root = validate_roots(input_root, output_root)
    sources = []
    for job in geometry_jobs(
        input_root, study=study, run_ids=run_ids, combinations=False
    ):
        manifest = _status(job, output_root)
        if manifest is None:
            raise ValueError(
                f"Missing or stale lobulus: {job.run_id}; run --step lobulus"
            )
        path = geometry_path(job.directory(output_root))
        source = job.sources[0]
        sources.append(
            RunSource(
                output_root,
                cast(Study, job.study),
                job.run_id,
                source.condition,
                source.condition_source,
                path,
                path.parent / "model.feb",
            )
        )
    return sources
