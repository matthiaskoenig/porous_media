"""Build an offline inspection report from already published SPT figures."""

from __future__ import annotations

import argparse
import csv
import json
import sys
import uuid
from collections import Counter
from dataclasses import dataclass, field
from datetime import UTC, datetime
from pathlib import Path
from typing import cast
from urllib.parse import quote

from jinja2 import Environment, FileSystemLoader, StrictUndefined

from research import RESULTS_DIR
from research.spt.artifacts import atomic_text, generation_path, validate_output

SECTIONS = {
    "analysis": "Quantitative analysis",
    "overviews": "Geometry overviews",
    "combinations": "Combined sector layouts",
    "geometry": "Geometry time courses",
    "supplements": "Supplemental geometry figures",
    "conditions": "Simulation conditions",
}
IMAGE_FORMATS = {".png", ".svg", ".pdf", ".gif", ".mp4"}


@dataclass
class Figure:
    """Describe one figure, grouping its exports into a single report card."""

    title: str
    section: str
    family: str
    geometry: str
    scale: str
    caption: str
    links: dict[str, str]
    patterns: list[str] = field(default_factory=list)
    doses: list[str] = field(default_factory=list)
    preview: str = ""
    manifest: str = ""
    provenance: str = ""
    identifier: str = ""
    keywords: str = ""
    order: int = 0


def _mapping(value: object) -> dict[str, object]:
    return cast(dict[str, object], value) if isinstance(value, dict) else {}


def _read_json(path: Path) -> dict[str, object]:
    return _mapping(json.loads(path.read_text(encoding="utf-8")))


def _local(path: Path, root: Path) -> bool:
    """Reject symbolic links anywhere below the chosen output root."""
    try:
        relative = path.relative_to(root)
    except ValueError:
        return False
    return not any(part in {"..", "."} for part in relative.parts) and not any(
        parent.is_symlink()
        for parent in (path, *path.parents)
        if parent != root and root in parent.parents
    )


def _href(path: Path, root: Path) -> str:
    return quote(path.relative_to(root).as_posix(), safe="/")


def _outputs(manifest: Path, root: Path) -> tuple[dict[str, object], list[Path]]:
    """Read only the active generation; allow relocation of the whole results tree.

    Output records retain their original absolute prefix after relocation. Only
    their suffix within the declared UUID generation is mapped locally. Source
    records are never followed, so the raw share can be offline.
    """
    if not _local(manifest, root):
        raise ValueError("Manifest is a symbolic link")
    value = _read_json(manifest)
    directory = generation_path(manifest.parent, value)
    name = directory.name
    if not _local(directory, root):
        raise ValueError("Missing or symbolic-link generation")
    records = value.get("outputs")
    if not isinstance(records, list) or not records:
        raise ValueError("Missing published outputs")
    outputs = []
    for record in records:
        record = _mapping(record)
        parts = Path(str(record.get("path", ""))).parts
        if parts.count(name) != 1 or ".." in parts:
            raise ValueError("Output outside its declared generation")
        path = directory.joinpath(*parts[parts.index(name) + 1 :])
        if not _local(path, root) or not path.is_file():
            raise ValueError(f"Missing or symbolic-link output: {path.name}")
        stat = path.stat()
        if stat.st_size != record.get("size") or stat.st_mtime_ns != record.get(
            "mtime_ns"
        ):
            raise ValueError(f"Changed output: {path.name}; rerun its processing step")
        outputs.append(path)
    return value, outputs


def _catalogue(root: Path) -> list[dict[str, str]]:
    path = root / "catalog.tsv"
    if not _local(path, root) or not path.is_file():
        return []
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _dose(row: dict[str, str]) -> str:
    value = row.get("boundary_flow_key", "")
    return str(int(float(value)) + 1) if value else ""


def _description(stem: str, spatial_times: list[object]) -> str:
    linear = stem.removesuffix("-symlog")
    titles = {
        "timecourses": "Matched-time courses",
        "endpoints": "Matched-time endpoints",
        "convergence-errors": "Differences from the finest available mesh",
        "initial-properties": "Initial cell properties",
        "final-state": "Final state overview",
        "meshes": "Distinct meshes",
        "timecourse": "Geometry time course",
    }
    if linear.startswith("spatial-"):
        token = linear.removeprefix("spatial-")
        if token.isdigit() and int(token) < len(spatial_times):
            return f"Cell values at {float(str(spatial_times[int(token)])):g} s"
        return f"Cell values {token}"
    return titles.get(linear, linear.replace("-", " ").replace("_", " "))


def _figure_order(stem: str) -> int:
    """Keep time courses and spatial snapshots in their natural reading order."""
    name = stem.removesuffix("-symlog")
    if name.startswith("spatial-") and name[8:].isdigit():
        return 10 + int(name[8:])
    return {
        "timecourses": 0,
        "final-state": 0,
        "meshes": 1,
        "endpoints": 100,
        "convergence-errors": 110,
        "initial-properties": 120,
    }.get(name, 0)


def _collect(
    root: Path, rows: list[dict[str, str]]
) -> tuple[list[Figure], list[dict[str, object]], list[str], set[str]]:
    figures: list[Figure] = []
    downloads: list[dict[str, object]] = []
    warnings: list[str] = []
    analysed: set[str] = set()
    conditions = {row["run_id"]: row for row in rows}
    patterns = (
        "analysis/*/manifest.json",
        "overviews/*/manifest.json",
        "figures/*/*/manifest.json",
        "*_combinations/xdmf/*/manifest.json",
    )
    bases = [
        root,
        *sorted(
            path for path in root.iterdir() if _local(path, root) and path.is_dir()
        ),
    ]
    manifests = sorted(
        {path for base in bases for pattern in patterns for path in base.glob(pattern)}
    )
    seen: set[Path] = set()
    for manifest in manifests:
        try:
            value, paths = _outputs(manifest, root)
            configuration = _mapping(value.get("configuration"))
            relative = manifest.relative_to(root)
            kind = next(
                part
                for part in relative.parts
                if part in {"analysis", "overviews", "figures", "xdmf"}
            )
            study = str(configuration.get("study", ""))
            if not study:
                study = (
                    manifest.parents[2].name
                    if kind == "xdmf"
                    else manifest.parent.name
                    if kind in {"analysis", "overviews"}
                    else manifest.parent.parent.name
                )
            family = "Convergence" if study.startswith("convergence") else "Simulations"
            geometry = (
                "Combined"
                if "combinations" in study
                else "Sixth"
                if "sixth" in study
                else "Lobulus"
            )
            section = {
                "analysis": "analysis",
                "overviews": "overviews",
                "xdmf": "combinations",
                "figures": "geometry",
            }[kind]
            supplemental = (
                relative.parts[0] not in {"analysis", "overviews", "figures"}
                and kind != "xdmf"
            )
            if supplemental:
                section = "supplements"
            runs = configuration.get("runs", [])
            runs = [str(run) for run in runs] if isinstance(runs, list) else []
            if kind == "overviews":
                sources = configuration.get("sources", [])
                if isinstance(sources, list):
                    runs = [str(source).split(" | ")[0] for source in sources]
            if kind == "analysis":
                analysed.update(runs)
            run = manifest.parent.name if kind in {"figures", "xdmf"} else ""
            selected = [conditions[key] for key in runs if key in conditions]
            if run in conditions:
                selected = [conditions[run]]
            if run.startswith("dose_") and run[5:].isdigit():
                selected = [row for row in rows if _dose(row) == str(int(run[5:]))]
            elif family == "Simulations" and geometry == "Combined":
                selected = [row for row in rows if row.get("pattern_name")]
            patterns_here = sorted(
                {row.get("pattern_name", "") for row in selected} - {""}
            )
            doses = sorted({_dose(row) for row in selected} - {""}, key=int)
            if run.startswith("dose_") and run[5:].isdigit():
                doses = [str(int(run[5:]))]
            layers = _mapping(configuration.get("options")).get("layers", [])
            layers = layers if isinstance(layers, list) else []
            symlog = any(_mapping(layer).get("scale") == "symlog" for layer in layers)
            keywords = " ".join(
                str(_mapping(layer).get("title", "")) for layer in layers
            )
            if kind == "analysis":
                keywords += (
                    " substrate product protein toxicant necrosis pressure volume"
                )
            provenance = next(
                (path for path in paths if path.name == "provenance.json"), None
            )
            spatial_times: list[object] = []
            if kind == "analysis" and provenance:
                analyses = _read_json(provenance).get("analyses", [])
                if isinstance(analyses, list) and analyses:
                    times = _mapping(_mapping(analyses[0]).get("configuration")).get(
                        "spatial_times", []
                    )
                    if isinstance(times, list):
                        spatial_times = times
            groups: dict[Path, dict[str, str]] = {}
            for path in paths:
                if path.suffix in IMAGE_FORMATS:
                    groups.setdefault(path.with_suffix(""), {})[
                        path.suffix[1:].upper()
                    ] = _href(path, root)
            for stem, links in sorted(groups.items()):
                if stem in seen:
                    continue
                seen.add(stem)
                scale = (
                    "Symlog" if symlog or stem.name.endswith("-symlog") else "Linear"
                )
                title = str(configuration.get("title", "")) or _description(
                    stem.name, spatial_times
                )
                if kind == "xdmf":
                    title = (
                        f"Dose {doses[0]}: six zonation patterns"
                        if doses
                        else "Six convergence resolutions"
                    )
                if (
                    kind == "analysis"
                    and stem.name.removesuffix("-symlog") == "endpoints"
                    and family == "Convergence"
                    and spatial_times
                ):
                    reference = (
                        5400.0
                        if 5400.0 in spatial_times
                        else float(str(spatial_times[len(spatial_times) // 2]))
                    )
                    title = f"Mesh resolution at {reference / 3600:g} h"
                caption = " · ".join(
                    filter(
                        None,
                        [
                            family,
                            geometry,
                            run,
                            *patterns_here,
                            *(f"dose {dose}" for dose in doses),
                        ],
                    )
                )
                if kind == "analysis":
                    caption = f"{family} · homogeneous lobules · {len(runs)} runs"
                if supplemental:
                    caption += f" · {relative.parts[0]}"
                figures.append(
                    Figure(
                        title,
                        section,
                        family,
                        geometry,
                        scale,
                        caption,
                        links,
                        patterns_here,
                        doses,
                        links.get("PNG", links.get("GIF", "")),
                        _href(manifest, root),
                        _href(provenance, root) if provenance else "",
                        keywords=" ".join(
                            [
                                keywords,
                                *runs,
                                *patterns_here,
                                *(f"dose {dose}" for dose in doses),
                            ]
                        ),
                        order=_figure_order(stem.name),
                    )
                )
            if kind == "analysis":
                downloads.append(
                    {
                        "title": f"{family} analysis",
                        "links": {
                            path.name: _href(path, root)
                            for path in paths
                            if path.suffix in {".tsv", ".nc", ".json"}
                        },
                    }
                )
        except (OSError, ValueError, StopIteration) as error:
            warnings.append(f"{manifest.relative_to(root)}: {error}")
    for name, title in (
        ("boundary_color", "Pattern and flow colour mapping"),
        ("boundary_flux", "Boundary flux conditions"),
    ):
        path = root / "figures" / f"{name}.png"
        if _local(path, root) and path.is_file():
            href = _href(path, root)
            figures.append(
                Figure(
                    title,
                    "conditions",
                    "Simulations",
                    "Conditions",
                    "Linear",
                    "Source generator conditions; flow units unspecified",
                    {"PNG": href},
                    preview=href,
                )
            )
    order = list(SECTIONS)
    figures.sort(
        key=lambda figure: (
            order.index(figure.section),
            figure.family != "Simulations",
            figure.geometry,
            figure.caption,
            figure.order,
            figure.title,
            figure.scale,
        )
    )
    for figure in figures:
        filename = (
            next(iter(figure.links.values())).rsplit("/", 1)[-1].rsplit(".", 1)[0]
        )
        key = f"{figure.manifest}/{filename}"
        figure.identifier = f"figure-{uuid.uuid5(uuid.NAMESPACE_URL, key).hex[:12]}"
    return figures, downloads, warnings, analysed


def build_report(output_root: Path) -> Path:
    """Publish report.html and index.html without reading raw or numerical data.

    All published figure generations below the output root are included, even
    when the most recent pipeline invocation selected only a subset. Incomplete
    or modified generations are listed as omissions. Images use relative links;
    move the entire results tree together to retain offline access.
    """
    root = output_root.resolve()
    if not root.is_dir():
        raise ValueError(f"Results directory does not exist: {root}")
    rows = _catalogue(root)
    figures, downloads, warnings, analysed = _collect(root, rows)
    pipeline_path = root / "pipeline-report.json"
    pipeline = (
        _read_json(pipeline_path)
        if _local(pipeline_path, root) and pipeline_path.is_file()
        else {}
    )
    counts = Counter(figure.section for figure in figures)
    resources = []
    for pattern in ("*catalog*.tsv", "*report*.json", "*validation*.json"):
        resources.extend(
            path
            for path in sorted(root.glob(pattern))
            if _local(path, root) and path.is_file()
        )
    downloads.append(
        {
            "title": "Catalogues and execution records",
            "links": {path.name: _href(path, root) for path in sorted(set(resources))},
        }
    )
    environment = Environment(
        loader=FileSystemLoader(Path(__file__).parent / "report_assets"),
        autoescape=True,
        undefined=StrictUndefined,
    )
    document = environment.get_template("report.html").render(
        figures=figures,
        sections=[
            {"key": key, "title": title, "count": counts[key]}
            for key, title in SECTIONS.items()
            if counts[key]
        ],
        generated=datetime.now(UTC).strftime("%Y-%m-%d %H:%M UTC"),
        analysed=len(analysed),
        exports=sum(len(figure.links) for figure in figures),
        pipeline=pipeline,
        warnings=warnings,
        rows=rows,
        truncated=[row for row in rows if row.get("complete") == "False"],
        downloads=downloads,
        filters={
            "family": sorted({figure.family for figure in figures}),
            "geometry": sorted({figure.geometry for figure in figures}),
            "scale": sorted({figure.scale for figure in figures}),
            "pattern": sorted(
                {pattern for figure in figures for pattern in figure.patterns}
            ),
            "dose": sorted(
                {dose for figure in figures for dose in figure.doses}, key=int
            ),
        },
    )
    # Render everything before replacing either entry point. Never follow links
    # supplied as destinations, including links to other research output files.
    for name in ("report.html", "index.html"):
        validate_output(root / name)
    for name in ("report.html", "index.html"):
        atomic_text(root / name, document)
    return root / "report.html"


def main(argv: list[str] | None = None) -> int:
    """Refresh the final report independently of processing and the raw share."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", type=Path, default=RESULTS_DIR / "spt")
    args = parser.parse_args(argv)
    try:
        print(build_report(args.output_root))
    except (OSError, ValueError) as error:
        print(f"SPT report: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
