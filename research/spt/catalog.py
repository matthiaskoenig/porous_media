"""Discover SPT runs and inspect XDMF metadata without loading numerical arrays."""

from __future__ import annotations

import math
import re
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import numpy as np

from research.spt.conditions import (
    SimulationCondition,
    conditions_from_source,
    generation_script,
    pattern_name,
)

type Study = Literal[
    "simulations_sixth",
    "convergence_sixth",
    "simulations_lobulus",
    "convergence_lobulus",
    "simulations_combinations",
    "convergence_combinations",
]
type CatalogValue = str | int | float | bool | None
STUDIES: tuple[Study, ...] = ("simulations_sixth", "convergence_sixth")
CONVERGENCE_MARKER = "convergence_"


def mesh_key(run_id: str) -> str:
    """Return the mesh resolution a convergence run carries in its name.

    A study prefixes its convergence runs with the variant they were generated
    for, so the resolution is what follows the last ``convergence_`` marker
    rather than the whole name after a fixed prefix.
    """
    index = run_id.rfind(CONVERGENCE_MARKER)
    return run_id if index < 0 else run_id[index + len(CONVERGENCE_MARKER) :]


@dataclass(frozen=True)
class RunSource:
    """Describe a run and the metadata that supplies its conditions."""

    root: Path
    study: Study
    run_id: str
    condition: SimulationCondition | None
    condition_source: str
    descriptor: Path | None = None
    model: Path | None = None

    @property
    def xdmf(self) -> Path:
        """Return the input descriptor."""
        return (
            self.descriptor or self.root / self.study / "xdmf" / f"{self.run_id}.xdmf"
        )

    @property
    def material(self) -> Path:
        """Return the corresponding model or simulation template.

        Every simulation of a study is rendered from one template, whose name
        carries the boundary conditions it was built for and therefore changes
        with the study design.
        """
        if self.model:
            return self.model
        folder = self.root / self.study / "material"
        if not self.study.startswith("simulations_"):
            return folder / f"{self.run_id}.feb"
        templates = sorted(folder.glob("*.feb.template"))
        if len(templates) > 1:
            raise ValueError(
                f"Ambiguous simulation templates: "
                f"{', '.join(path.name for path in templates)}"
            )
        return (
            templates[0]
            if templates
            else folder / "lobule_zonation_pattern.feb.template"
        )


@dataclass
class ParsedRun:
    """Hold descriptor metadata and scalar field references for one run."""

    source: RunSource
    row: dict[str, CatalogValue]
    times: np.ndarray
    grids: list[ET.Element]
    references: dict[str, Path]
    sources: tuple[Path, ...]
    fields: dict[str, tuple[str, ...]]

    def count(self, center: str) -> int:
        """Return the validated number of cells or points."""
        value = self.row["num_cells" if center == "Cell" else "num_points"]
        if not isinstance(value, int):
            raise ValueError("Invalid mesh count")
        return value


def discover_runs(
    root: Path,
    studies: tuple[Study, ...] = STUDIES,
    run_ids: tuple[str, ...] = (),
) -> list[RunSource]:
    """List known study directories once and reject empty or unknown selections."""
    root = root.expanduser().resolve()
    conditions, origin = (
        conditions_from_source(root) if "simulations_sixth" in studies else ([], "")
    )
    by_id = {condition.sim_id: condition for condition in conditions}
    runs = []
    for study in studies:
        folder = root / study / "xdmf"
        if not folder.is_dir():
            raise FileNotFoundError(f"Missing study directory: {folder}")
        paths = sorted(folder.glob("*.xdmf"))
        if not paths:
            raise ValueError(f"No XDMF runs found in {folder}")
        for path in paths:
            if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*", path.stem):
                raise ValueError(f"Invalid run ID: {path.stem!r}")
            if run_ids and path.stem not in run_ids:
                continue
            if study == "simulations_sixth" and path.stem not in by_id:
                raise ValueError(f"No generation conditions for {path.stem}")
            runs.append(RunSource(root, study, path.stem, by_id.get(path.stem), origin))
    unknown = set(run_ids) - {run.run_id for run in runs}
    if unknown:
        raise ValueError(f"Unknown runs: {', '.join(sorted(unknown))}")
    return runs


def _parse_xml(path: Path) -> ET.ElementTree[ET.Element]:
    try:
        return ET.parse(path)
    except ET.ParseError as error:
        raise ValueError(f"Invalid XML in {path}: {error}") from error


def _parse_model(path: Path) -> ET.ElementTree[ET.Element]:
    """Parse a FEBio model, tolerating a comment XML rejects but FEBio accepts.

    One archived template delimits a commented-out load with ``<!---`` and
    ``--->``. XML forbids a comment whose content ends in a dash, so the file
    does not parse although FEBio ran it. Comments carry no model semantics,
    so they are removed and the model is parsed again; a file that still fails
    keeps its original error.
    """
    try:
        return _parse_xml(path)
    except ValueError as error:
        stripped = re.sub(rb"<!--.*?--+>", b"", path.read_bytes(), flags=re.DOTALL)
        try:
            return ET.ElementTree(ET.fromstring(stripped))
        except ET.ParseError:
            raise error from None


def expected_duration(material: Path) -> float | None:
    """Read one explicit FEBio control duration, leaving ambiguous models unknown."""
    if not material.is_file():
        return None
    controls = _parse_model(material).findall(".//Control")
    if len(controls) != 1:
        return None
    steps = controls[0].findtext("time_steps")
    size = controls[0].findtext("step_size")
    if steps is None or size is None:
        return None
    duration = float(steps) * float(size)
    if not math.isfinite(duration) or duration <= 0:
        raise ValueError(f"Invalid model duration: {material}")
    return duration


def hdf_reference(item: ET.Element) -> tuple[str, str]:
    """Require a simple relative HDF reference rather than an executable input."""
    if item.get("Format") != "HDF" or len(item):
        raise ValueError("SPT processing requires simple HDF-backed DataItems")
    filename, separator, dataset = (item.text or "").strip().partition(":")
    if not separator or not dataset.startswith("/") or not filename:
        raise ValueError("Invalid XDMF HDF reference")
    if Path(filename).is_absolute() or "\\" in filename:
        raise ValueError("HDF references must be relative to the XDMF descriptor")
    return filename, dataset


def hyperslab(item: ET.Element) -> tuple[ET.Element, tuple[slice, ...]]:
    """Split a HyperSlab into its HDF source and the region it selects.

    Stacked series keep one dataset per attribute and give every step a
    HyperSlab into it, so the selection is what separates the steps.
    """
    children = list(item)
    if len(children) != 2:
        raise ValueError("A HyperSlab needs a selection and a source DataItem")
    selection, source = children
    if selection.get("Format") != "XML" or len(selection):
        raise ValueError("HyperSlab selections must be inline XML")
    values = [int(value) for value in (selection.text or "").split()]
    if not values or len(values) % 3:
        raise ValueError("A HyperSlab selection needs start, stride and count")
    rank = len(values) // 3
    start, stride, count = values[:rank], values[rank : 2 * rank], values[2 * rank :]
    if any(begin < 0 for begin in start) or any(step <= 0 for step in stride + count):
        raise ValueError("HyperSlab selections must be positive and start inside")
    region = tuple(
        slice(begin, begin + size * step, step)
        for begin, step, size in zip(start, stride, count, strict=True)
    )
    return source, region


def hdf_source(item: ET.Element) -> tuple[str, str, tuple[slice, ...] | None]:
    """Resolve a data item, which may select one step of a stacked attribute."""
    if item.get("ItemType") != "HyperSlab":
        filename, dataset = hdf_reference(item)
        return filename, dataset, None
    source, region = hyperslab(item)
    filename, dataset = hdf_reference(source)
    return filename, dataset, region


def hdf_sources(tree: ET.ElementTree[ET.Element]) -> list[ET.Element]:
    """List the data items naming HDF datasets, descending into HyperSlabs."""
    sources: list[ET.Element] = []
    nested: set[int] = set()
    for item in tree.iter("DataItem"):
        if id(item) in nested:
            continue
        if item.get("ItemType") == "HyperSlab":
            source, _ = hyperslab(item)
            nested.update(id(child) for child in item)
            sources.append(source)
        else:
            sources.append(item)
    return sources


def inspect_run(source: RunSource) -> ParsedRun:
    """Validate time axes and scalar schemas using descriptor and model metadata."""
    tree = _parse_xml(source.xdmf)
    collections = tree.findall(".//Grid[@CollectionType='Temporal']")
    if len(collections) != 1:
        raise ValueError("Expected one XDMF temporal collection")
    grids = collections[0].findall("Grid")
    if not grids:
        raise ValueError("Empty XDMF time series")
    time_values = []
    for grid in grids:
        time = grid.find("Time")
        if time is None:
            raise ValueError("Missing time value")
        time_values.append(float(time.get("Value", "nan")))
    times = np.array(time_values)
    if not np.isfinite(times).all() or np.any(np.diff(times) <= 0):
        raise ValueError("Time values must be finite and strictly increasing")
    topology = tree.find(".//Topology")
    geometry = tree.find(".//Geometry/DataItem")
    if topology is None or geometry is None:
        raise ValueError("Missing mesh topology or geometry")
    counts = {
        "Cell": int(topology.attrib["NumberOfElements"]),
        "Node": int(geometry.attrib["Dimensions"].split()[0]),
    }
    if any(count <= 0 for count in counts.values()):
        raise ValueError("SPT meshes must contain points and cells")
    fields: dict[str, tuple[str, ...]] = {}
    for grid in grids:
        current: dict[str, list[str]] = {"Cell": [], "Node": []}
        for attr in grid.findall("Attribute"):
            center = attr.get("Center", "Node")
            item = attr.find("DataItem")
            if center not in current or item is None:
                raise ValueError("Unsupported XDMF attribute")
            dimensions = tuple(
                int(value) for value in item.attrib["Dimensions"].split()
            )
            if attr.get("AttributeType", "Scalar") == "Scalar" and dimensions not in {
                (counts[center],),
                (counts[center], 1),
            }:
                raise ValueError(
                    f"Scalar dimensions differ from the mesh: {attr.attrib['Name']}"
                )
            if dimensions in {(counts[center],), (counts[center], 1)}:
                current[center].append(attr.attrib["Name"])
        schema = {key: tuple(sorted(names)) for key, names in current.items()}
        if any(len(set(names)) != len(names) for names in current.values()):
            raise ValueError("Duplicate scalar field names")
        if fields and schema != fields:
            raise ValueError("Scalar field schema changes between time steps")
        fields = schema
    references: dict[str, Path] = {}
    for item in hdf_sources(tree):
        filename, _ = hdf_reference(item)
        if filename not in references:
            path = (source.xdmf.parent / filename).resolve()
            if not path.is_relative_to(source.root) or not path.is_file():
                raise ValueError(
                    f"HDF sidecar missing or outside input root: {filename}"
                )
            references[filename] = path
    duration = expected_duration(source.material)
    condition = source.condition
    row: dict[str, CatalogValue] = {
        "study": source.study,
        "run_id": source.run_id,
        "xdmf": str(source.xdmf),
        "num_cells": counts["Cell"],
        "num_points": counts["Node"],
        "num_steps": len(times),
        "tstart_s": float(times[0]),
        "tend_s": float(times[-1]),
        "expected_duration_s": duration,
        "complete": None
        if duration is None
        else bool(times[-1] - times[0] >= duration - 1e-6),
        "cell_fields": "|".join(fields["Cell"]),
        "point_fields": "|".join(fields["Node"]),
        "resolution_key": mesh_key(source.run_id) if condition is None else None,
        "pattern_key": None if condition is None else condition.pattern_key,
        "pattern_name": None
        if condition is None
        else pattern_name(condition.pattern_key),
        "boundary_flow_key": None if condition is None else condition.boundary_flow_key,
        "boundary_flow": None if condition is None else condition.boundary_flow,
        "boundary_flow_unit": "unspecified in source" if condition else None,
        "condition_source": source.condition_source
        if condition
        else str(source.material),
    }
    sources = {source.xdmf, source.material, *references.values()}
    script = generation_script(source.root)
    if condition:
        sources.add(script)
    return ParsedRun(
        source, row, times, grids, references, tuple(sorted(sources)), fields
    )
