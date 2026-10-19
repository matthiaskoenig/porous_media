"""Process supplied SPT XDMF/NetCDF results into resumable local analysis caches."""

import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict
from multiprocessing import get_context
from pathlib import Path
from typing import cast

import pandas as pd

from porous_media.data.xdmf_tools import xdmfs_from_directory
from research.spt.cache import (
    LegacyPolicy,
    ProcessingResult,
    Profile,
    cache_status,
    process_run,
    validate_roots,
)
from research.spt.catalog import CatalogValue, discover_runs, inspect_run
from research.spt.cli import (
    atomic_json,
    selected_studies,
    selection_parser,
    write_catalog,
)
from research.spt.conditions import conditions_from_source


def process_febio_simulations(input_dir: Path, xdmf_dir: Path) -> dict[Path, Path]:
    """Convert earlier VTK exports explicitly into a separate output directory."""
    validate_roots(input_dir, xdmf_dir)
    return xdmfs_from_directory(input_dir=input_dir, xdmf_dir=xdmf_dir, overwrite=False)


def process_studies(
    input_root: Path,
    output_root: Path,
    *,
    study: str = "both",
    run_ids: tuple[str, ...] = (),
    profile: Profile = "analysis",
    legacy: LegacyPolicy = "sample",
    workers: int = 2,
    force: bool = False,
    dry_run: bool = False,
) -> list[ProcessingResult]:
    """Discover, process, and report independent runs with bounded worker memory.

    A dry run reads metadata only and never creates output files or directories.
    One worker executes serially. Each worker owns one run and its HDF handles.
    """
    if workers < 1:
        raise ValueError("workers must be positive")
    source_root, target = validate_roots(input_root, output_root)
    runs = discover_runs(source_root, selected_studies(study), run_ids)
    if not run_ids and study in {"both", "simulations"}:
        conditions, _ = conditions_from_source(source_root)
        missing = {item.sim_id for item in conditions} - {
            source.run_id for source in runs if source.study == "simulations_sixth"
        }
        if missing:
            raise ValueError(
                f"Missing expected simulation runs: {', '.join(sorted(missing))}. "
                "Use --runs to process an explicit subset."
            )
    results = []
    if dry_run:
        for source in runs:
            try:
                status, manifest = cache_status(source, target, profile, legacy)
                row = (
                    cast(dict[str, CatalogValue], manifest["catalog"])
                    if manifest is not None and status == "current"
                    else inspect_run(source).row
                )
            except (OSError, ValueError, KeyError) as error:
                results.append(
                    ProcessingResult(
                        source.study, source.run_id, "failed", str(error), 0.0, {}
                    )
                )
                continue
            results.append(
                ProcessingResult(
                    source.study,
                    source.run_id,
                    status,
                    "rebuild requested" if force else status,
                    0.0,
                    row,
                )
            )
        return results
    pending = []
    for source in runs:
        status, manifest = cache_status(source, target, profile, legacy)
        if status == "current" and not force and manifest is not None:
            results.append(
                ProcessingResult(
                    source.study,
                    source.run_id,
                    "cached",
                    "cached",
                    0.0,
                    cast(dict[str, CatalogValue], manifest["catalog"]),
                )
            )
        else:
            pending.append(source)
    if workers == 1 or len(pending) <= 1:
        for source in pending:
            result = process_run(source, target, profile, legacy, force)
            results.append(result)
            print(
                f"{result.study}/{result.run_id}: {result.status} ({result.elapsed_s:.2f}s) {result.message}",
                flush=True,
            )
    else:
        with ProcessPoolExecutor(
            max_workers=min(workers, len(pending)), mp_context=get_context("spawn")
        ) as executor:
            futures = {
                executor.submit(
                    process_run, source, target, profile, legacy, force
                ): source
                for source in pending
            }
            for future in as_completed(futures):
                try:
                    result = future.result()
                except Exception as error:
                    source = futures[future]
                    result = ProcessingResult(
                        source.study,
                        source.run_id,
                        "failed",
                        f"Worker failed: {error}",
                        0.0,
                        {},
                    )
                results.append(result)
                print(
                    f"{result.study}/{result.run_id}: {result.status} ({result.elapsed_s:.2f}s) {result.message}",
                    flush=True,
                )
    results.sort(key=lambda result: (result.study, result.run_id))
    atomic_json(
        target / "processing-report.json", [asdict(result) for result in results]
    )
    write_catalog(
        pd.DataFrame([result.row for result in results if result.row]),
        target,
        name="processing-catalog",
    )
    return results


def main(argv: list[str] | None = None) -> int:
    """Run the local processing CLI and return nonzero if any run fails."""
    parser = selection_parser(__doc__ or "Process SPT data")
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--profile", choices=("analysis", "full"), default="analysis")
    parser.add_argument(
        "--legacy-cache",
        choices=("sample", "full", "ignore"),
        default="sample",
        help="Validate first/middle/last frames, verify every frame, or rebuild from HDF",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Rebuild selected runs from primary HDF data",
    )
    parser.add_argument(
        "--dry-run",
        "--status",
        action="store_true",
        help="Inspect metadata and cache state without writing",
    )
    args = parser.parse_args(argv)
    try:
        results = process_studies(
            args.input_root,
            args.output_root,
            study=args.study,
            run_ids=tuple(args.runs),
            profile=args.profile,
            legacy=args.legacy_cache,
            workers=args.workers,
            force=args.force,
            dry_run=args.dry_run,
        )
    except (OSError, ValueError) as error:
        print(f"SPT processing: {error}", file=sys.stderr)
        return 1
    if args.dry_run:
        for result in results:
            print(
                f"{result.study}/{result.run_id}: {result.message}; end={result.row.get('tend_s')}s; complete={result.row.get('complete')}"
            )
    counts = {
        state: sum(result.status == state for result in results)
        for state in sorted({result.status for result in results})
    }
    print(f"{len(results)} runs: {counts}")
    truncated = [
        result.run_id for result in results if result.row.get("complete") is False
    ]
    if truncated:
        print(
            f"Runs ending before model duration: {', '.join(truncated)}",
            file=sys.stderr,
        )
    return int(any(result.status == "failed" for result in results))


if __name__ == "__main__":
    raise SystemExit(main())
