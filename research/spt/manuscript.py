"""Create annotated SPT manuscript Figures 4-7 and a convergence supplement."""

from __future__ import annotations

import argparse
import sys
import tempfile
from dataclasses import asdict
from pathlib import Path

import matplotlib as mpl
import numpy as np

from research.spt import results_spt_dir
from research.spt.analysis_plots import PlotOptions, _curves
from research.spt.artifacts import (
    atomic_json,
    atomic_output,
    atomic_text,
    current_manifest,
    fingerprint,
    generation_path,
    publish,
    run_lock,
    validate_output,
)
from research.spt.fields import SUMMARY_FIELDS
from research.spt.layout import (
    category_root,
    exports_directory,
    log_path,
    manuscript_directory,
)
from research.spt.manuscript_data import (
    CONVERGENCE_DISPLAY_MESHES,
    CONVERGENCE_TIME,
    SIMULATION_TIMES,
    ManuscriptData,
    pattern_label,
)
from research.spt.manuscript_geometry import (
    render_example,
    render_necrosis,
    render_patterns,
    save_manuscript_figure,
)
from research.spt.manuscript_patterns import PATTERN_ICON_DESCRIPTION
from research.spt.manuscript_timecourses import (
    CONVERGENCE_COLORS,
    annotate_timecourses,
    render_convergence,
)
from research.spt.plot_data import PlotSource
from research.spt.plot_render import RenderOptions
from research.spt.plot_view import GEOMETRY_VIEW
from research.spt.report import build_report

MANUSCRIPT_VERSION = 4
PAPER = "https://doi.org/10.1101/2024.03.26.586870"
FIGURES = ("4", "5", "6", "7", "S1")
TITLES = {
    "4": "Figure 4 | SPT metabolism over time",
    "5": "Figure 5 | Zonation pattern comparison",
    "6": "Figure 6 | Mean and SD time courses",
    "7": "Figure 7 | Necrosis by pattern and substrate inflow",
    "S1": "Supplement | Convergence meshes and time courses",
}


def export_figure(generation: Path, root: Path, stem: str) -> None:
    """Expose stable manuscript filenames without duplicating generation storage."""
    directory = exports_directory(root)
    validate_output(directory / stem)
    directory.mkdir(parents=True, exist_ok=True)
    for name in (
        f"{stem}.svg",
        f"{stem}.pdf",
        f"{stem}.png",
        "caption.md",
        "provenance.json",
    ):
        source = generation / name
        target = directory / (
            f"{stem}.{name}" if name in {"caption.md", "provenance.json"} else name
        )
        validate_output(target)
        if not source.is_file():
            target.unlink(missing_ok=True)
        elif not target.exists() or not target.samefile(source):
            with atomic_output(target) as temporary:
                temporary.unlink()
                temporary.hardlink_to(source)


def simulation_cohort(data: ManuscriptData) -> list[PlotSource]:
    """Require all six patterns and eight doses, identified by recorded conditions."""
    sources = [
        data.source("simulations_lobulus", path.parent.name)
        for path in sorted(
            (category_root(data.root, "analysis-data") / "simulations_lobulus").glob(
                "*/manifest.json"
            )
        )
    ]
    conditions = [
        (source.run.condition.pattern_key, source.run.condition.boundary_flow_key)
        for source in sources
        if source.run.condition
    ]
    if len(sources) != 48 or set(conditions) != {
        (p, d) for p in range(6) for d in range(8)
    }:
        raise ValueError(
            "Manuscript figures require six patterns and eight doses; run the complete SPT pipeline"
        )
    return sorted(
        sources,
        key=lambda source: (
            source.run.condition.pattern_key if source.run.condition else -1,
            source.run.condition.boundary_flow_key if source.run.condition else -1,
        ),
    )


def caption(
    number: str,
    sources: list[PlotSource],
    dose: int,
    comparison_time: float = CONVERGENCE_TIME,
) -> str:
    """Describe the actual run selection, timing, normalization and statistics."""
    if number == "4":
        return (
            f"Figure 4. SPT metabolism for the linear increasing zonation pattern, "
            f"D{dose} ({sources[0].run.run_id}). Rows show extracellular "
            "substrate, extracellular product, protein, toxic compound and necrosis. "
            "Columns show the actual 0, 6, 12, 18 and 24 h snapshots from left to "
            "right. Continuous cell "
            "values are interpolated between native frames; discrete necrosis holds "
            "the preceding frame. One shared vertical colorbar beside each field row uses scalar limits "
            "from full histories of all 48 runs. Maps show the lobule outline "
            "without internal mesh edges."
        )
    if number == "5":
        return (
            f"Figure 5. SPT metabolism at D{dose} for A: constant, B: linear "
            "increase, C: sharp perivenous, D: linear decrease, E: sharp periportal "
            "and F: random zonation. Each block has three rows for protein, toxic "
            "compound and necrosis. Columns show 0, 6, 12, 18 and 24 h from left to "
            "right. All blocks share full-cohort scalar limits and a single bottom "
            "set of three field colorbars. Maps show lobule outlines without "
            "internal mesh edges. Runs: "
            + ", ".join(source.run.run_id for source in sources)
            + "."
        )
    if number == "6":
        return (
            "Figure 6. Volume-weighted mean and population SD over homogeneous "
            "lobules for all 48 simulations. Rows A-F use the same pattern order as "
            "Figure 5. Eight shades per pattern encode increasing doses D1-D8. "
            "Necrosis is expressed as percent tissue volume. Means and SD reuse "
            "the existing matched-time analysis, which interpolates cell values and "
            "volumes before reduction. SD describes spatial heterogeneity, not "
            "uncertainty across independent simulations."
        )
    if number == "7":
        return (
            "Figure 7. Necrosis at 24 h by protein zonation pattern (rows A-F) and "
            "increasing substrate inflow (D1-D8, columns). The shared binary legend "
            "on the right labels black as necrosis and white as no necrosis. "
            "Cell states retain the fixed 0-1 scale. Maps show lobule outlines "
            "without internal mesh edges."
        )
    common_end = min(source.end for source in sources)
    return (
        "Convergence supplement. A: source-sector meshes with their triangles, "
        "ordered coarse to fine, "
        "plus the mixed lobule with sectors 1-6 clockwise from upper left "
        "(37250, 25000, 18750, 12500, 09375, 06250). "
        "Mesh boundaries and native time courses share one distinct color per "
        "resolution. The 50000 mesh is excluded from this figure. "
        f"B: mixed-sixth substrate, product, protein, toxic compound and necrosis "
        f"at the common output time {comparison_time:g} s ({comparison_time / 3600:g} h). "
        "Canonical field palettes and full-history limits shared across the six "
        "meshes permit direct sector comparisons. Only sector outlines are shown. "
        "Continuous values use native-frame interpolation; necrosis holds its "
        "preceding native state. C: native volume-weighted lobule time courses "
        "with linear axes. Crosses "
        "mark the actual ends of truncated histories, without extrapolation. The "
        f"dotted line marks the shared interval ending at {common_end / 3600:.4f} h. "
        "Mesh identifiers are the archived file labels; cell counts refer to sixths."
    )


def create_manuscript_figures(
    data_root: Path,
    *,
    figures: tuple[str, ...] | None = None,
    dose: int = 6,
    formats: tuple[str, ...] = ("svg", "pdf", "png"),
    dpi: int = 300,
    comparison_time: float = CONVERGENCE_TIME,
    force: bool = False,
) -> dict[str, Path]:
    """Publish figures atomically from local caches and reuse unchanged exports.

    The original manuscript's representative dose is index 6. Snapshot labels
    are generated from actual seconds, independently of the overview selections.
    Existing analyses, geometry caches and general-purpose figures are not changed.
    """
    options = RenderOptions(formats=formats, dpi=dpi)
    options.validate()
    if figures is None:
        figures = FIGURES
    if not figures or len(set(figures)) != len(figures) or set(figures) - set(FIGURES):
        raise ValueError("Select unique figure identifiers: 4, 5, 6, 7, S1")
    if not 1 <= dose <= 8:
        raise ValueError("Dose index must be between 1 and 8")
    if not np.isfinite(comparison_time) or comparison_time < 0:
        raise ValueError("Comparison time must be finite and nonnegative")
    data = ManuscriptData(data_root)
    simulations = simulation_cohort(data) if set(figures) - {"S1"} else []
    chosen = [source for source in simulations if source.dose == dose - 1]
    ranges = data.ranges(simulations) if simulations else {}
    output = {}
    for number in figures:
        if number == "S1":
            sources = [
                data.source("convergence_lobulus", f"convergence_{mesh}")
                for mesh in CONVERGENCE_DISPLAY_MESHES
            ]
            if (
                not max(source.start for source in sources)
                <= comparison_time
                <= min(source.end for source in sources)
            ):
                raise ValueError(
                    "Convergence comparison time is outside the common native history"
                )
        elif number == "4":
            sources = [
                source
                for source in chosen
                if source.run.condition and source.run.condition.pattern_key == 1
            ]
        else:
            sources = chosen if number == "5" else simulations
        mesh_sources = (
            []
            if number == "6"
            else [(source.run.study, source.run.run_id) for source in sources]
        )
        if number == "S1":
            mesh_sources.append(("convergence_combinations", "convergence_mixed"))
        inputs = data.dependencies(
            simulations if number in {"4", "5"} else sources, mesh_sources
        )
        configuration: dict[str, object] = {
            "version": MANUSCRIPT_VERSION,
            "figure": number,
            "title": TITLES[number],
            "study": sources[0].run.study,
            "runs": [source.run.run_id for source in sources],
            "paper": PAPER,
            "dose_index": dose if number in {"4", "5"} else None,
            "snapshot_times_s": list(SIMULATION_TIMES)
            if number in {"4", "5"}
            else [86400.0]
            if number == "7"
            else [],
            "formats": list(formats),
            "dpi": dpi,
            "crop": "visible content with a one-point safety border",
            "scalar_ranges": {key: list(value) for key, value in ranges.items()}
            if number in {"4", "5"}
            else {},
            "mesh_rendering": "raster at export DPI; editable vector annotations and curves",
        }
        stem = f"figure-{int(number):02}" if number != "S1" else "figure-s1-convergence"
        rendered_stem = stem
        if number != "6":
            configuration["geometry_view"] = GEOMETRY_VIEW
        if number in {"4", "5", "6", "7"}:
            configuration.update(
                manuscript_layout_version=3,
                dose_labels=[f"D{index}" for index in range(1, 9)],
                pattern_labels=[pattern_label(index) for index in range(6)],
                pattern_icons=PATTERN_ICON_DESCRIPTION,
                colorbars="five shared vertical field bars"
                if number == "4"
                else "three shared horizontal field bars"
                if number == "5"
                else "one shared binary necrosis legend"
                if number == "7"
                else "none",
            )
        if number in {"4", "5"}:
            configuration.update(
                time_grid_version=2,
                time_axis="columns",
                map_fields=list(
                    SUMMARY_FIELDS if number == "4" else SUMMARY_FIELDS[2:]
                ),
            )
        if number == "S1":
            configuration.update(
                convergence_readouts_version=1,
                convergence_layout_version=3,
                timecourse_scale="linear",
                snapshot_times_s=[comparison_time],
                mixed_sector_meshes=list(CONVERGENCE_DISPLAY_MESHES),
                mixed_sector_direction="clockwise from upper left",
                mesh_colors=CONVERGENCE_COLORS,
                scalar_ranges={
                    field: list(limits)
                    for field, limits in data.ranges(sources).items()
                },
            )
        directory = manuscript_directory(data.root, stem)
        with run_lock(directory):
            cached = current_manifest(directory, configuration, data.root)
            if cached and not force:
                output[number] = generation_path(directory, cached)
                export_figure(output[number], data.root, stem)
                print(f"Figure {number}: cached", flush=True)
                continue
            records = [fingerprint(path) for path in inputs]
            with tempfile.TemporaryDirectory(
                prefix=".partial-", dir=directory
            ) as temporary:
                staging = Path(temporary)
                with mpl.rc_context({"svg.fonttype": "none", "pdf.fonttype": 42}):
                    if number == "4":
                        render_example(
                            data,
                            sources[0],
                            ranges,
                            staging,
                            options,
                        )
                    elif number == "5":
                        render_patterns(data, sources, ranges, staging, options)
                    elif number == "7":
                        render_necrosis(data, sources, staging, options)
                    elif number == "S1":
                        render_convergence(
                            data,
                            sources,
                            staging,
                            options,
                            comparison_time=comparison_time,
                        )
                    else:
                        datasets = [data.summary(source) for source in sources]
                        try:
                            _curves(
                                sources,
                                datasets,
                                staging,
                                stem,
                                PlotOptions(formats=formats, dpi=dpi),
                                SUMMARY_FIELDS,
                                decorate=annotate_timecourses,
                                export=save_manuscript_figure,
                            )
                        finally:
                            for dataset in datasets:
                                dataset.close()
                if stem != rendered_stem and number != "6":
                    for extension in formats:
                        (staging / f"{rendered_stem}.{extension}").rename(
                            staging / f"{stem}.{extension}"
                        )
                description = caption(number, sources, dose, comparison_time)
                if number != "6":
                    description += (
                        " Geometry uses an orthographic view tilted 20 degrees from "
                        "above (azimuth -60 degrees). Native thickness is retained "
                        "without exaggeration; caps and exterior walls share their "
                        "owning cell's unmodified scalar color."
                    )
                if number in {"4", "5", "6", "7"}:
                    description += " " + PATTERN_ICON_DESCRIPTION
                atomic_text(staging / "caption.md", description + "\n")
                atomic_json(
                    staging / "provenance.json",
                    {
                        **configuration,
                        "caption": description,
                        "conditions": [
                            asdict(source.run.condition)
                            if source.run.condition
                            else None
                            for source in sources
                        ],
                        "native_end_s": {
                            source.run.run_id: source.end for source in sources
                        },
                        "raw_hdf_reads": 0,
                    },
                )
                output[number] = publish(
                    directory, staging, configuration, inputs, records
                )
            export_figure(output[number], data.root, stem)
            print(f"Figure {number}: rendered", flush=True)
    return output


def main(argv: list[str] | None = None) -> int:
    """Create selected manuscript figures from a completed local SPT pipeline."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, default=results_spt_dir)
    parser.add_argument("--figures", choices=FIGURES, nargs="+")
    parser.add_argument("--dose", type=int, default=6)
    parser.add_argument(
        "--formats",
        choices=("svg", "pdf", "png"),
        nargs="+",
        default=["svg", "pdf", "png"],
    )
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument(
        "--comparison-time",
        type=float,
        default=CONVERGENCE_TIME,
        help="Common convergence map time in seconds (default: 5400, or 1.5 h)",
    )
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args(argv)
    try:
        outputs = create_manuscript_figures(
            args.data_root,
            figures=tuple(args.figures) if args.figures is not None else None,
            dose=args.dose,
            formats=tuple(args.formats),
            dpi=args.dpi,
            comparison_time=args.comparison_time,
            force=args.force,
        )
        atomic_json(
            log_path(args.data_root, "manuscript.json"),
            {
                "status": "passed",
                "figures": {key: str(path) for key, path in outputs.items()},
                "raw_hdf_reads": 0,
            },
        )
        print(build_report(args.data_root))
    except (OSError, ValueError, KeyError) as error:
        print(f"Manuscript figures failed: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
