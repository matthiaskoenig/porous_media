"""Run resumable quantitative SPT analysis on homogeneous lobules by default."""

from __future__ import annotations

import html
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict
from multiprocessing import get_context
from pathlib import Path
from typing import cast

import pandas as pd

from research import RESULTS_DIR
from research.spt.analysis_data import AnalysisOptions, prepare_analysis, timelines
from research.spt.analysis_plots import (
    PlotOptions,
    cached_analysis_figures,
    render_analysis,
)
from research.spt.cache import ProcessingResult, validate_roots
from research.spt.cli import atomic_json, selection_parser, write_catalog
from research.spt.plot_data import PlotSource, plot_sources


def _prepare(
    source: PlotSource,
    data_root: Path,
    output_root: Path,
    options: AnalysisOptions,
    times: tuple[float, ...],
    spatial: tuple[float, ...],
    force: bool,
) -> tuple[dict[str, object] | None, str]:
    try:
        manifest = prepare_analysis(
            source, data_root, output_root, options, times, spatial, force=force
        )
        if manifest is None:
            raise RuntimeError(
                "Analysis inputs changed; retry after processing completes"
            )
        return manifest, ""
    except Exception as error:
        return None, f"{type(error).__name__}: {error}"


def process_analysis(
    input_root: Path,
    data_root: Path,
    output_root: Path,
    *,
    study: str = "both",
    geometry: str = "lobulus",
    run_ids: tuple[str, ...] = (),
    options: AnalysisOptions | None = None,
    plots: PlotOptions | None = None,
    workers: int = 2,
    force: bool = False,
    dry_run: bool = False,
) -> list[ProcessingResult]:
    """Reduce each local cache once and share its results between study figures."""
    started = time.perf_counter()
    input_root, output_root = validate_roots(input_root, output_root)
    _, data_root = validate_roots(input_root, data_root)
    options = options or AnalysisOptions()
    plots = plots or PlotOptions()
    options.validate()
    plots.validate()
    if workers < 1 or geometry not in {"sixth", "lobulus"}:
        raise ValueError(
            "Use positive workers and homogeneous sixth or lobulus geometry"
        )
    sources = plot_sources(
        input_root, data_root, study=study, geometry=geometry, run_ids=run_ids
    )
    if any(output_root.is_relative_to(source.run.xdmf.parent) for source in sources):
        raise ValueError("Analysis outputs must not modify source history directories")
    windows = timelines(sources, options)
    prepared: dict[str, dict[str, object]] = {}
    results = []
    pending = []
    for source in sources:
        try:
            times, spatial = windows[source.family]
            cached = prepare_analysis(
                source, data_root, output_root, options, times, spatial, dry_run=True
            )
            if cached and not force:
                prepared[source.run.run_id] = cached
                results.append(
                    ProcessingResult(
                        source.run.study,
                        source.run.run_id,
                        "cached",
                        "analysis data current",
                        0.0,
                        {},
                    )
                )
            elif dry_run:
                results.append(
                    ProcessingResult(
                        source.run.study,
                        source.run.run_id,
                        "pending",
                        "analysis data required",
                        0.0,
                        {},
                    )
                )
            else:
                pending.append(source)
        except (OSError, ValueError, KeyError) as error:
            results.append(
                ProcessingResult(
                    source.run.study, source.run.run_id, "failed", str(error), 0.0, {}
                )
            )

    def record(
        source: PlotSource, outcome: tuple[dict[str, object] | None, str]
    ) -> None:
        manifest, error = outcome
        if manifest is not None:
            prepared[source.run.run_id] = manifest
        status = "failed" if error else "processed"
        results.append(
            ProcessingResult(
                source.run.study,
                source.run.run_id,
                status,
                error or "analysis data prepared",
                0.0,
                {},
            )
        )
        print(
            f"analysis {source.run.study}/{source.run.run_id}: {status} {error}",
            flush=True,
        )

    if pending and workers > 1:
        with ProcessPoolExecutor(
            max_workers=min(workers, len(pending)), mp_context=get_context("spawn")
        ) as executor:
            futures = {
                executor.submit(
                    _prepare,
                    source,
                    data_root,
                    output_root,
                    options,
                    *windows[source.family],
                    force,
                ): source
                for source in pending
            }
            for future in as_completed(futures):
                source = futures[future]
                try:
                    outcome = future.result()
                except Exception as error:
                    outcome = None, f"Worker failed: {error}"
                record(source, outcome)
    else:
        for source in pending:
            record(
                source,
                _prepare(
                    source,
                    data_root,
                    output_root,
                    options,
                    *windows[source.family],
                    force,
                ),
            )
    results.sort(key=lambda result: (result.study, result.run_id))
    if dry_run:
        for family in sorted({source.family for source in sources}):
            selected = [source for source in sources if source.family == family]
            if any(
                source.run.run_id not in prepared for source in selected
            ) or not cached_analysis_figures(selected, output_root, plots):
                results.append(
                    ProcessingResult(
                        selected[0].run.study,
                        "figures",
                        "pending",
                        "study figures or serialized tables require generation",
                        0.0,
                        {},
                    )
                )
        return results
    rows = []
    for family in sorted({source.family for source in sources}):
        selected = [source for source in sources if source.family == family]
        if any(source.run.run_id not in prepared for source in selected):
            continue
        try:
            generation = render_analysis(
                selected,
                [prepared[source.run.run_id] for source in selected],
                output_root,
                plots,
                force=force,
            )
            row = {
                "study": selected[0].run.study,
                "runs": len(selected),
                "gallery": str(generation / "index.html"),
                "netcdf": str(generation / "analysis.nc"),
                "native": str(generation / "native-timecourses.tsv"),
                "matched": str(generation / "matched-timecourses.tsv"),
                "profiles": str(generation / "spatial-profiles.tsv"),
            }
            rows.append(row)
            print(f"figures {row['study']}: ready", flush=True)
        except Exception as error:
            results.append(
                ProcessingResult(
                    selected[0].run.study,
                    "figures",
                    "failed",
                    f"{type(error).__name__}: {error}",
                    0.0,
                    {},
                )
            )
    directory = output_root / "analysis"
    if directory.resolve() != directory.absolute():
        raise ValueError("Analysis index must not traverse symbolic links")
    directory.mkdir(parents=True, exist_ok=True)
    links = "".join(
        f'<li><a href="{html.escape(Path(str(row["gallery"])).relative_to(directory).as_posix(), quote=True)}">{html.escape(str(row["study"]))}</a> ({row["runs"]} runs)</li>'
        for row in rows
    )
    temporary = directory / ".index.partial"
    if temporary.is_symlink() or (directory / "index.html").is_symlink():
        raise ValueError("Analysis index must not follow symbolic links")
    temporary.write_text(
        f'<!doctype html><html lang="en"><meta charset="utf-8"><title>SPT quantitative analysis</title><h1>SPT quantitative analysis</h1><ul>{links}</ul><p>Time courses retain native histories. Matched comparisons interpolate cell values and volumes before reduction. SD describes spatial variation.</p></html>',
        encoding="utf-8",
    )
    temporary.replace(directory / "index.html")
    write_catalog(pd.DataFrame(rows), output_root, name=f"analysis-catalog-{study}")
    atomic_json(
        output_root / f"analysis-report-{study}.json",
        {
            "elapsed_s": time.perf_counter() - started,
            "results": [asdict(result) for result in results],
            "outputs": rows,
            "native_frames_read": sum(
                cast(
                    int,
                    cast(dict[str, object], manifest["metadata"])["native_frames_read"],
                )
                for key, manifest in prepared.items()
                if any(
                    result.run_id == key and result.status == "processed"
                    for result in results
                )
            ),
            "windows": windows,
            "geometry": geometry,
            "data_root": str(data_root),
        },
    )
    return results


def main(argv: list[str] | None = None, *, study: str = "both") -> int:
    """Analyze local caches and return nonzero for missing inputs or failed outputs."""
    parser = selection_parser(__doc__ or "Analyze SPT")
    parser.set_defaults(study=study)
    parser.add_argument("--data-root", type=Path, default=RESULTS_DIR / "spt")
    parser.add_argument("--geometry", choices=("sixth", "lobulus"), default="lobulus")
    parser.add_argument("--fields", nargs="+", default=list(AnalysisOptions().fields))
    parser.add_argument("--time-points", type=int, default=51)
    parser.add_argument("--bins", type=int, default=20)
    parser.add_argument("--spatial-times", nargs="+", type=float, default=[])
    parser.add_argument(
        "--formats",
        nargs="+",
        choices=("svg", "pdf", "png"),
        default=["svg", "pdf", "png"],
    )
    parser.add_argument("--dpi", type=int, default=110)
    parser.add_argument(
        "--spread",
        action="store_true",
        help="Show volume-weighted population SD, not confidence intervals",
    )
    parser.add_argument(
        "--linear-only",
        action="store_true",
        help="Omit supplemental convergence symlog figures",
    )
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--status", "--dry-run", action="store_true")
    args = parser.parse_args(argv)
    try:
        results = process_analysis(
            args.input_root,
            args.data_root,
            args.output_root,
            study=args.study,
            geometry=args.geometry,
            run_ids=tuple(args.runs),
            workers=args.workers,
            options=AnalysisOptions(
                tuple(args.fields),
                args.time_points,
                args.bins,
                tuple(args.spatial_times),
            ),
            plots=PlotOptions(
                tuple(dict.fromkeys(args.formats)),
                args.dpi,
                args.spread,
                not args.linear_only,
            ),
            force=args.force,
            dry_run=args.status,
        )
    except (OSError, ValueError, KeyError) as error:
        print(f"SPT analysis: {error}", file=sys.stderr)
        return 1
    for result in results:
        if args.status or result.status == "failed":
            print(f"{result.study}/{result.run_id}: {result.status} {result.message}")
    counts = {
        state: sum(result.status == state for result in results)
        for state in sorted({result.status for result in results})
    }
    print(f"{len(results)} analysis runs: {counts}")
    return int(any(result.status == "failed" for result in results))


if __name__ == "__main__":
    raise SystemExit(main())
