"""Read publication inputs from validated local SPT caches without raw HDF reads."""

from __future__ import annotations

import json
from collections.abc import Sequence
from pathlib import Path
from typing import cast

import numpy as np
import xarray as xr

from research.spt.analysis_data import (
    ANALYSIS_VERSION,
    CellSampler,
    input_paths,
    load_analysis,
)
from research.spt.artifacts import current_manifest, generation_path
from research.spt.cache import _current_generation, open_processed_cells
from research.spt.catalog import RunSource, Study
from research.spt.conditions import SimulationCondition
from research.spt.fields import SUMMARY_FIELDS
from research.spt.geometry import CONVERGENCE_MESHES
from research.spt.plot_data import SAMPLE_VERSION, PlotSource
from research.spt.plot_render import Samples

SIMULATION_TIMES = (0.0, 21600.0, 43200.0, 64800.0, 86400.0)
CONVERGENCE_TIME = 5400.0
PATTERN_NAMES = (
    "constant",
    "linear increase",
    "sharp perivenous",
    "linear decrease",
    "sharp periportal",
    "random",
)


class PublicationData:
    """Share small summaries, geometry arrays and selected cell frames between figures."""

    def __init__(self, root: Path) -> None:
        """Use an existing pipeline results tree; no network input root is required."""
        self.root = root.expanduser().resolve()
        self._sources: dict[tuple[str, str], PlotSource] = {}
        self._meshes: dict[str, Samples] = {}
        self._snapshots: dict[tuple[str, tuple[float, ...]], Samples] = {}

    def source(self, study: str, run_id: str) -> PlotSource:
        """Validate a serialized homogeneous-lobule analysis and its dependencies."""
        key = (study, run_id)
        if key in self._sources:
            return self._sources[key]
        if study not in {"simulations_lobulus", "convergence_lobulus"}:
            raise ValueError(
                f"Publication analysis requires homogeneous lobules: {study}"
            )
        directory = self.root / "analysis-data" / study / run_id
        manifest = json.loads((directory / "manifest.json").read_text())
        configuration, metadata = manifest["configuration"], manifest["metadata"]
        if (
            configuration.get("version") != ANALYSIS_VERSION
            or current_manifest(directory, configuration, self.root) is None
        ):
            raise ValueError(f"Analysis is stale; rerun the SPT pipeline: {run_id}")
        descriptor = Path(configuration["source"])
        if not descriptor.resolve().is_relative_to(self.root):
            raise ValueError("Publication inputs must use local derived geometries")
        condition = configuration.get("condition")
        selected = (
            SimulationCondition(
                str(condition["sim_id"]),
                int(condition["pattern_key"]),
                int(condition["boundary_flow_key"]),
                float(condition["boundary_flow"]),
            )
            if condition is not None
            else None
        )
        if selected is not None and (
            selected.sim_id != run_id
            or selected.pattern_key not in range(6)
            or selected.boundary_flow_key not in range(8)
        ):
            raise ValueError(f"Invalid publication run condition: {run_id}")
        source = PlotSource(
            RunSource(
                self.root,
                cast(Study, study),
                run_id,
                selected,
                "validated local analysis",
                descriptor,
                descriptor.parent / "model.feb",
            ),
            "lobulus",
            float(metadata["source_start_s"]),
            float(metadata["source_end_s"]),
            str(metadata["title"]),
            dose=selected.boundary_flow_key if selected else None,
        )
        self._sources[key] = source
        return source

    def mesh_manifest(self, study: str, run_id: str) -> tuple[Path, dict[str, object]]:
        """Check cached plotting geometry without opening a mesh or numerical array."""
        directory = self.root / "plot-data" / study / run_id
        manifest = json.loads((directory / "manifest.json").read_text())
        configuration = manifest["configuration"]
        if (
            configuration.get("version") != SAMPLE_VERSION
            or current_manifest(directory, configuration, self.root) is None
        ):
            raise ValueError(
                f"Geometry samples are stale; rerun geometry plotting: {run_id}"
            )
        return directory, manifest

    def mesh(self, study: str, run_id: str) -> Samples:
        """Load each distinct projected mesh once, keeping native cell identities."""
        directory, manifest = self.mesh_manifest(study, run_id)
        metadata = cast(dict[str, object], manifest["metadata"])
        mesh_hash = str(metadata["mesh_hash"])
        if mesh_hash not in self._meshes:
            samples = Samples.load(generation_path(directory, manifest) / "samples.npz")
            # The geometry cache can be shared across runs; its old scalar frames cannot.
            self._meshes[mesh_hash] = Samples(
                samples.faces,
                samples.cell_ids,
                np.empty((0, 0, 0)),
                samples.sectors,
                np.empty(0),
            )
        return self._meshes[mesh_hash]

    def snapshots(self, source: PlotSource, times: tuple[float, ...]) -> Samples:
        """Interpolate local continuous cells and hold discrete necrosis at true times."""
        if (
            not times
            or not np.isfinite(times).all()
            or np.any(np.diff(times) <= 0)
            or min(times) < source.start
            or max(times) > source.end
        ):
            raise ValueError(
                f"Snapshot times outside native history: {source.run.run_id}"
            )
        key = (source.run.run_id, times)
        if key not in self._snapshots:
            mesh = self.mesh(source.run.study, source.run.run_id)
            with open_processed_cells(source.run, self.root) as cells:
                sampler = CellSampler(cells, SUMMARY_FIELDS)
                frames = [
                    np.stack([frame[name] for name in SUMMARY_FIELDS])
                    for frame in (sampler.sample(time) for time in times)
                ]
            values = np.stack(frames)
            if values.shape[-1] != len(mesh.sectors):
                raise ValueError("Geometry and processed cell identities differ")
            self._snapshots[key] = Samples(
                mesh.faces,
                mesh.cell_ids,
                values,
                mesh.sectors,
                np.array(times),
            )
        return self._snapshots[key]

    def summary(self, source: PlotSource, *, native: bool = False) -> xr.Dataset:
        """Load only compact volume-weighted summaries, retaining actual native ends."""
        return load_analysis(source, self.root, "native" if native else "aligned")

    def mixed_snapshot(self, time: float) -> Samples:
        """Assemble mixed sixths from native local frames at one common time.

        Composition preserves cell order within each sector. Validate its source
        order and cell counts before transferring scalar values from sector zero
        of each homogeneous lobule. Never interpolate already sampled plot frames.
        """
        study, run_id = "convergence_combinations", "convergence_mixed"
        _, manifest = self.mesh_manifest(study, run_id)
        configuration = cast(dict[str, object], manifest["configuration"])
        descriptor = Path(str(configuration["descriptor"]))
        if not descriptor.resolve().is_relative_to(self.root):
            raise ValueError("Mixed comparison geometry must be local")
        composition = json.loads(
            (descriptor.parent.parent / "manifest.json").read_text()
        )
        recorded = [Path(path).stem for path in composition["configuration"]["sources"]]
        expected = [f"convergence_{mesh}" for mesh in CONVERGENCE_MESHES]
        if recorded != expected:
            raise ValueError("Mixed convergence sector source order differs")
        mixed = self.mesh(study, run_id)
        if set(np.unique(mixed.sectors)) != set(range(6)):
            raise ValueError("Mixed convergence requires exactly six sectors")
        values = np.empty((1, len(SUMMARY_FIELDS), len(mixed.sectors)))
        for sector, name in enumerate(expected):
            source = self.source("convergence_lobulus", name)
            samples = self.snapshots(source, (time,))
            selected, target = samples.sectors == 0, mixed.sectors == sector
            if np.count_nonzero(selected) != np.count_nonzero(target):
                raise ValueError(f"Mixed convergence cell identities differ: {name}")
            values[0][:, target] = samples.values[0][:, selected]
        return Samples(
            mixed.faces, mixed.cell_ids, values, mixed.sectors, np.array([time])
        )

    def ranges(self, sources: list[PlotSource]) -> dict[str, tuple[float, float]]:
        """Reuse full-history cell limits so changing snapshots does not change colors."""
        ranges = {field: (float("inf"), float("-inf")) for field in SUMMARY_FIELDS}
        for source in sources:
            directory = _current_generation(source.run, self.root)
            values = json.loads((directory / "limits.json").read_text())["Cell"]
            for field, (low, high) in ranges.items():
                ranges[field] = (
                    min(low, float(values[field]["min"])),
                    max(high, float(values[field]["max"])),
                )
        if not np.isfinite(list(ranges.values())).all():
            raise ValueError("Publication scalar limits must be finite")
        return ranges

    def dependencies(
        self, sources: list[PlotSource], meshes: Sequence[tuple[str, str]]
    ) -> tuple[Path, ...]:
        """Fingerprint stable manifest pointers after validating their local outputs."""
        paths = {
            source.directory(self.root, "analysis-data") / "manifest.json"
            for source in sources
        }
        for source in sources:
            paths.update(input_paths(source, self.root))
            paths.add(_current_generation(source.run, self.root) / "limits.json")
        for study, run_id in meshes:
            directory, _ = self.mesh_manifest(study, run_id)
            paths.add(directory / "manifest.json")
        return tuple(sorted(paths))
