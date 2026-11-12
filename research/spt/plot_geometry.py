"""Create fast vector figures for SPT sixths, lobules and six-sector comparisons."""

from __future__ import annotations

import json
import sys
import time
from contextlib import suppress
from dataclasses import asdict, dataclass, replace
from functools import partial
from pathlib import Path
from typing import cast

import numpy as np
import pandas as pd

from research.spt import results_spt_dir
from research.spt.artifacts import atomic_json
from research.spt.cache import ProcessingResult, _current_generation, validate_roots
from research.spt.cli import selection_parser, write_catalog
from research.spt.execution import execute_jobs, print_result
from research.spt.fields import LAYERS, SPT_LAYERS, Layer
from research.spt.layout import log_path
from research.spt.plot_data import (
    PlotSource,
    plot_sources,
    prepare_samples,
    sample_times,
)
from research.spt.plot_overviews import overview_figures
from research.spt.plot_render import RenderOptions, render_source


@dataclass
class Prepared:
    """Report selective sampling independently from figure generation."""

    source: PlotSource
    manifest: dict[str, object] | None
    error: str = ""
    created: bool = False


def _prepare(
    source: PlotSource,
    output: Path,
    fields: tuple[str, ...],
    times: dict[str, tuple[float, ...]],
    force: bool,
) -> Prepared:
    try:
        return Prepared(
            source,
            prepare_samples(source, output, fields, times[source.family], force=force),
            created=True,
        )
    except Exception as error:
        return Prepared(source, None, f"{type(error).__name__}: {error}")


def _render(
    item: Prepared,
    output: Path,
    options: dict[str, RenderOptions],
    ranges: dict[str, dict[str, tuple[float, float]]],
    force: bool,
) -> ProcessingResult:
    return render_source(
        item.source,
        cast(dict[str, object], item.manifest),
        output,
        options[item.source.family],
        ranges[item.source.family],
        force=force,
    )


def shared_ranges(
    prepared: list[Prepared],
    data_root: Path | None = None,
) -> dict[str, dict[str, tuple[float, float]]]:
    """Share full-history source limits, falling back to sampled synthetic data."""
    result: dict[str, dict[str, tuple[float, float]]] = {}
    for item in prepared:
        if item.manifest is None:
            continue
        metadata = cast(dict[str, object], item.manifest["metadata"])
        ranges = cast(dict[str, list[float]], metadata["ranges"])
        limits_path = item.source.run.xdmf.with_name(
            f"{item.source.run.xdmf.stem}_limits.json"
        )
        full = {}
        if data_root is not None:
            try:
                cache = _current_generation(item.source.run, data_root)
                cached_limits = json.loads(
                    (cache / "limits.json").read_text(encoding="utf-8")
                )["Cell"]
                full = {
                    name: [values["min"], values["max"]]
                    for name, values in cached_limits.items()
                }
            except (OSError, ValueError, KeyError):
                pass
        if set(ranges) - full.keys() and limits_path.is_file():
            with suppress(OSError, ValueError, KeyError):
                full = {
                    **json.loads(limits_path.read_text(encoding="utf-8"))["limits"],
                    **full,
                }
        ranges = {name: full.get(name, values) for name, values in ranges.items()}
        selected = result.setdefault(item.source.family, {})
        for name, values in ranges.items():
            low, high = selected.get(name, (float("inf"), float("-inf")))
            selected[name] = (min(low, values[0]), max(high, values[1]))
    return result


def process_plots(
    input_root: Path,
    data_root: Path,
    output_root: Path,
    *,
    study: str = "both",
    geometry: str = "all",
    run_ids: tuple[str, ...] = (),
    snapshots: int = 5,
    times: tuple[float, ...] = (),
    options: RenderOptions | None = None,
    workers: int = 2,
    overviews: bool = True,
    per_history: bool = True,
    video_frames: int = 0,
    force: bool = False,
    dry_run: bool = False,
    provided_sources: list[PlotSource] | None = None,
) -> list[ProcessingResult]:
    """Sample requested data once, render selected outputs, and resume completed work."""
    started = time.perf_counter()
    input_root, output_root = validate_roots(input_root, output_root)
    _, data_root = validate_roots(input_root, data_root)
    if workers < 1:
        raise ValueError("workers must be positive")
    options = options or RenderOptions()
    if not per_history and not overviews:
        raise ValueError("Select per-history figures or overviews")
    sources = provided_sources or plot_sources(
        input_root,
        data_root,
        study=study,
        geometry=geometry,
        run_ids=run_ids,
        plot_root=output_root,
    )
    if any(output_root.is_relative_to(source.run.xdmf.parent) for source in sources):
        raise ValueError("Plot outputs must not modify a source history directory")
    timelines = sample_times(
        sources, snapshots, times, video_frames, manuscript_snapshots=True
    )
    if not per_history and (video_frames or options.panels):
        raise ValueError("Videos and standalone panels require per-history rendering")
    family_options = {
        family: replace(
            options,
            snapshots=values[0],
            video_times=values[1],
        )
        for family, values in timelines.items()
    }
    for selected in family_options.values():
        selected.validate()
    fields = tuple(layer.sid for layer in options.layers)
    all_times = {
        family: tuple(
            float(value) for value in np.union1d(chosen.snapshots, chosen.video_times)
        )
        for family, chosen in family_options.items()
    }
    prepared = []
    pending = []
    for source in sources:
        manifest = prepare_samples(
            source, output_root, fields, all_times[source.family], dry_run=True
        )
        if manifest is not None and not force:
            prepared.append(Prepared(source, manifest))
        elif dry_run:
            prepared.append(Prepared(source, None))
        else:
            pending.append(source)
    for source, item in execute_jobs(
        pending,
        partial(
            _prepare, output=output_root, fields=fields, times=all_times, force=force
        ),
        workers=workers,
        failed=lambda source, error: Prepared(source, None, f"Worker failed: {error}"),
    ):
        prepared.append(item)
        print(
            f"sample {source.run.study}/{source.run.run_id}: {'failed: ' + item.error if item.error else 'prepared'}",
            flush=True,
        )
    prepared.sort(key=lambda item: (item.source.run.study, item.source.run.run_id))
    ranges = shared_ranges(prepared, data_root)
    results = []
    render_jobs = []
    for item in prepared:
        source = item.source
        if item.manifest is None:
            results.append(
                ProcessingResult(
                    source.run.study,
                    source.run.run_id,
                    "failed" if item.error else "pending",
                    item.error or "sampling required",
                    0.0,
                    {},
                )
            )
            continue
        if not per_history:
            results.append(
                ProcessingResult(
                    source.run.study,
                    source.run.run_id,
                    "sampled",
                    "overview input ready",
                    0.0,
                    {},
                )
            )
            continue
        selected = family_options[source.family]
        status = render_source(
            source,
            item.manifest,
            output_root,
            selected,
            ranges[source.family],
            force=force,
            dry_run=True,
        )
        if dry_run:
            results.append(status)
        elif status.status == "cached" and not force:
            results.append(
                render_source(
                    source, item.manifest, output_root, selected, ranges[source.family]
                )
            )
        else:
            render_jobs.append(item)
    for _, result in execute_jobs(
        render_jobs,
        partial(
            _render,
            output=output_root,
            options=family_options,
            ranges=ranges,
            force=force,
        ),
        workers=workers,
        failed=lambda item, error: ProcessingResult.failed(
            item.source.run.study, item.source.run.run_id, error
        ),
    ):
        results.append(result)
        print_result(result, prefix="figure ")
    results.sort(key=lambda result: (result.study, result.run_id))
    if dry_run:
        return results
    overview_rows: list[dict[str, str]] = []
    if overviews:
        for name in sorted({item.source.run.study for item in prepared}):
            selected_items = [
                item
                for item in prepared
                if item.source.run.study == name and item.manifest is not None
            ]
            if not selected_items:
                continue
            try:
                overview_rows.extend(
                    overview_figures(
                        [item.source for item in selected_items],
                        [
                            cast(dict[str, object], item.manifest)
                            for item in selected_items
                        ],
                        output_root,
                        family_options[selected_items[0].source.family],
                        ranges[selected_items[0].source.family],
                        force=force,
                    )
                )
            except Exception as error:
                results.append(
                    ProcessingResult(name, "overviews", "failed", str(error), 0.0, {})
                )
    rows: list[dict[str, object]] = [
        dict(result.row) for result in results if result.row
    ]
    rows.extend(dict(row) for row in overview_rows)
    write_catalog(pd.DataFrame(rows), log_path(output_root, f"plots-{study}.tsv"))
    reads = sum(
        cast(int, cast(dict[str, object], item.manifest["metadata"])["hdf_field_reads"])
        for item in prepared
        if item.created and item.manifest is not None
    )
    atomic_json(
        log_path(output_root, f"plots-{study}.json"),
        {
            "elapsed_s": time.perf_counter() - started,
            "sampling": {
                "prepared": sum(item.created for item in prepared),
                "cached": sum(not item.created for item in prepared),
                "hdf_field_reads": reads,
            },
            "results": [asdict(result) for result in results],
            "overviews": overview_rows,
            "ranges": ranges,
        },
    )
    return results


def main(argv: list[str] | None = None, *, study: str = "both") -> int:
    """Render selected SPT geometry figures and return nonzero on any failure."""
    parser = selection_parser(__doc__ or "Plot SPT geometry")
    parser.set_defaults(study=study)
    parser.add_argument(
        "--data-root",
        type=Path,
        default=results_spt_dir,
        help="Root containing published lobulus and combination histories",
    )
    parser.add_argument(
        "--geometry", choices=("all", "sixth", "lobulus", "combined"), default="all"
    )
    parser.add_argument("--fields", nargs="+", default=[layer.sid for layer in LAYERS])
    parser.add_argument(
        "--scale",
        choices=("linear", "symlog"),
        default="linear",
        help="Continuous scalar color scale; necrosis retains a linear scale",
    )
    parser.add_argument(
        "--formats",
        nargs="+",
        choices=("svg", "pdf", "png"),
        default=["png"],
    )
    parser.add_argument("--snapshots", type=int, default=5)
    parser.add_argument(
        "--times",
        nargs="+",
        type=float,
        default=[],
        help="Explicit snapshot times in seconds, within every selected source history",
    )
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument(
        "--dpi",
        type=int,
        default=150,
        help="PNG preview and video resolution; SVG/PDF stay vector",
    )
    parser.add_argument(
        "--video-frames",
        type=int,
        default=0,
        help="Opt-in video frame count; no intermediate PNG files are written",
    )
    parser.add_argument("--fps", type=int, default=12)
    parser.add_argument("--gif", action="store_true")
    parser.add_argument(
        "--panels",
        action="store_true",
        help="Also save standalone requested field/time panels",
    )
    parser.add_argument("--no-overviews", action="store_true")
    parser.add_argument("--only-overviews", action="store_true")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--dry-run", "--status", action="store_true")
    args = parser.parse_args(argv)
    known = {layer.sid: layer for layer in SPT_LAYERS}
    layers = tuple(known.get(name, Layer(name, name)) for name in args.fields)
    if args.scale == "symlog":
        layers = tuple(
            replace(layer, scale="symlog", title=f"{layer.title} (symlog)")
            if layer.sid != "rr_necrosis"
            else layer
            for layer in layers
        )
    options = RenderOptions(
        layers,
        tuple(dict.fromkeys(args.formats)),
        dpi=args.dpi,
        fps=args.fps,
        panels=args.panels,
        gif=args.gif,
    )
    try:
        results = process_plots(
            args.input_root,
            args.data_root,
            args.output_root,
            study=args.study,
            geometry=args.geometry,
            run_ids=tuple(args.runs),
            snapshots=args.snapshots,
            times=tuple(args.times),
            options=options,
            workers=args.workers,
            overviews=not args.no_overviews,
            per_history=not args.only_overviews,
            video_frames=args.video_frames,
            force=args.force,
            dry_run=args.dry_run,
        )
    except (OSError, ValueError, KeyError) as error:
        print(f"SPT geometry plots: {error}", file=sys.stderr)
        return 1
    for result in results:
        if args.dry_run or result.status == "failed":
            print(f"{result.study}/{result.run_id}: {result.status} {result.message}")
    counts = {
        state: sum(result.status == state for result in results)
        for state in sorted({result.status for result in results})
    }
    print(f"{len(results)} histories: {counts}")
    return int(any(result.status == "failed" for result in results))


if __name__ == "__main__":
    raise SystemExit(main())
