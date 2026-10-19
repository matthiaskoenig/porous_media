"""Run the complete SPT catalogue, reconstruction, analysis and figure workflow."""

from __future__ import annotations

import sys
import time
from collections.abc import Callable
from dataclasses import asdict
from pathlib import Path
from typing import TypedDict

from research.spt.analysis import process_analysis
from research.spt.analysis_data import AnalysisOptions
from research.spt.analysis_plots import PlotOptions
from research.spt.cache import ProcessingResult, validate_roots
from research.spt.cli import atomic_json, selection_parser
from research.spt.data_processing import process_studies
from research.spt.geometry import process_geometries
from research.spt.information import main as information_main
from research.spt.plot_geometry import process_plots
from research.spt.plot_render import RenderOptions
from research.spt.report import build_report


class PipelineReport(TypedDict):
    """Describe the selected workflow and each completed or pending stage."""

    status: str
    elapsed_s: float
    input_root: str
    output_root: str
    study: str
    runs: list[str]
    stages: list[dict[str, object]]
    truncated_runs: list[str]
    scope: str


def run_pipeline(
    input_root: Path,
    output_root: Path,
    *,
    study: str = "both",
    run_ids: tuple[str, ...] = (),
    workers: int = 2,
    formats: tuple[str, ...] = ("svg", "pdf", "png"),
    analysis_options: AnalysisOptions | None = None,
    plot_options: PlotOptions | None = None,
    geometry_figures: bool = True,
    geometry_snapshots: int = 5,
    dry_run: bool = False,
) -> PipelineReport:
    """Resume every stage, stop on failure, and record a reviewable execution report.

    Existing generations are reused only when their source and output fingerprints
    match. This processes existing simulation results; it does not run FEBio.
    Explicit run subsets omit mixed geometries, which require complete studies.
    """
    started = time.perf_counter()
    input_root, output_root = validate_roots(input_root, output_root)
    analysis_options = analysis_options or AnalysisOptions()
    plot_options = plot_options or PlotOptions(formats=formats)
    analysis_options.validate()
    plot_options.validate()
    RenderOptions(formats=formats).validate()
    if workers < 1 or geometry_snapshots < 1:
        raise ValueError("workers and geometry_snapshots must be positive")
    stages: list[dict[str, object]] = []
    truncated: set[str] = set()
    failed = False

    def stage(name: str, action: Callable[[], list[ProcessingResult] | int]) -> bool:
        stage_started = time.perf_counter()
        print(f"SPT stage: {name}", flush=True)
        try:
            result = action()
            if not isinstance(result, int):
                truncated.update(
                    item.run_id for item in result if item.row.get("complete") is False
                )
            bad = (
                bool(result)
                if isinstance(result, int)
                else any(item.status == "failed" for item in result)
            )
            pending = (
                False
                if isinstance(result, int)
                else any(
                    item.status not in {"cached", "current", "processed", "rendered"}
                    for item in result
                )
            )
            stages.append(
                {
                    "name": name,
                    "status": "failed"
                    if bad
                    else "pending"
                    if dry_run and pending
                    else "passed",
                    "elapsed_s": time.perf_counter() - stage_started,
                    "results": result
                    if isinstance(result, int)
                    else [asdict(item) for item in result],
                }
            )
            return not bad
        except Exception as error:
            stages.append(
                {
                    "name": name,
                    "status": "failed",
                    "elapsed_s": time.perf_counter() - stage_started,
                    "error": f"{type(error).__name__}: {error}",
                }
            )
            return False

    args = [
        "--input-root",
        str(input_root),
        "--output-root",
        str(output_root),
        "--study",
        study,
    ]
    if run_ids:
        args += ["--runs", *run_ids]
    args += ["--dry-run"] if dry_run else ["--plots"]
    if not stage("catalogue", lambda: information_main(args)) or not stage(
        "geometry",
        lambda: process_geometries(
            input_root,
            output_root,
            study=study,
            run_ids=run_ids,
            combinations=not run_ids,
            workers=workers,
            dry_run=dry_run,
        ),
    ):
        failed = True
    elif dry_run and stages[-1]["status"] == "pending":
        stages.extend(
            {
                "name": name,
                "status": "pending",
                "reason": "requires completed geometries",
            }
            for name in ("numeric-cache", "analysis", "geometry-figures")
            if geometry_figures or name != "geometry-figures"
        )
    elif not stage(
        "numeric-cache",
        lambda: process_studies(
            input_root,
            output_root,
            study=study,
            run_ids=run_ids,
            geometry="lobulus",
            workers=workers,
            dry_run=dry_run,
        ),
    ):
        failed = True
    elif dry_run and stages[-1]["status"] == "pending":
        stages.append(
            {
                "name": "analysis",
                "status": "pending",
                "reason": "requires completed numerical caches",
            }
        )
        if geometry_figures:
            stages.append(
                {
                    "name": "geometry-figures",
                    "status": "pending",
                    "reason": "not inspected before quantitative analysis is ready",
                }
            )
    elif not stage(
        "analysis",
        lambda: process_analysis(
            input_root,
            output_root,
            output_root,
            study=study,
            run_ids=run_ids,
            options=analysis_options,
            plots=plot_options,
            workers=workers,
            dry_run=dry_run,
        ),
    ) or (
        geometry_figures
        and not stage(
            "geometry-figures",
            lambda: process_plots(
                input_root,
                output_root,
                output_root,
                study=study,
                run_ids=run_ids,
                options=RenderOptions(formats=formats),
                snapshots=geometry_snapshots,
                workers=workers,
                dry_run=dry_run,
            ),
        )
    ):
        failed = True
    report: PipelineReport = {
        "status": "failed"
        if failed
        else "pending"
        if any(item["status"] == "pending" for item in stages)
        else "passed",
        "elapsed_s": time.perf_counter() - started,
        "input_root": str(input_root),
        "output_root": str(output_root),
        "study": study,
        "runs": list(run_ids),
        "stages": stages,
        "truncated_runs": sorted(truncated),
        "scope": "existing solver results; homogeneous lobulus analysis; mixed sectors for visualization only",
    }
    if not dry_run:
        atomic_json(output_root / "pipeline-report.json", report)
        build_report(output_root)
    return report


def main(argv: list[str] | None = None) -> int:
    """Run all stages or inspect their status without changing inputs."""
    parser = selection_parser(__doc__ or "Run SPT")
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument(
        "--formats",
        nargs="+",
        choices=("svg", "pdf", "png"),
        default=["svg", "pdf", "png"],
    )
    parser.add_argument("--time-points", type=int, default=51)
    parser.add_argument("--bins", type=int, default=20)
    parser.add_argument("--spatial-times", nargs="+", type=float, default=[])
    parser.add_argument("--no-geometry-figures", action="store_true")
    parser.add_argument("--geometry-snapshots", type=int, default=5)
    parser.add_argument("--status", "--dry-run", action="store_true")
    args = parser.parse_args(argv)
    try:
        report = run_pipeline(
            args.input_root,
            args.output_root,
            study=args.study,
            run_ids=tuple(args.runs),
            workers=args.workers,
            formats=tuple(dict.fromkeys(args.formats)),
            analysis_options=AnalysisOptions(
                time_points=args.time_points,
                bins=args.bins,
                spatial_times=tuple(args.spatial_times),
            ),
            geometry_figures=not args.no_geometry_figures,
            geometry_snapshots=args.geometry_snapshots,
            dry_run=args.status,
        )
    except (OSError, ValueError, KeyError) as error:
        print(f"SPT pipeline: {error}", file=sys.stderr)
        return 1
    for item in report["stages"]:
        print(f"{item['name']}: {item['status']}")
        if "error" in item:
            print(item["error"], file=sys.stderr)
    print(f"SPT pipeline: {report['status']}")
    if report["truncated_runs"]:
        print(
            f"Source histories ending before model duration: {', '.join(report['truncated_runs'])}"
        )
    return int(report["status"] == "failed")


if __name__ == "__main__":
    raise SystemExit(main())
