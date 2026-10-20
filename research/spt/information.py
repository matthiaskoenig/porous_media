"""Catalogue SPT simulation/convergence runs and export their conditions locally."""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

from research.spt.conditions import (
    SimulationCondition,
    simulation_conditions,
)
from research.spt.conditions import boundary_flows as boundary_flows
from research.spt.conditions import pattern_idx2name as pattern_idx2name
from research.spt.plot_style import PATTERN_COLORMAPS, simulation_color

pattern_colormaps = PATTERN_COLORMAPS


def simulation_conditions_df(
    conditions: list[SimulationCondition] | None = None,
) -> pd.DataFrame:
    """Create a typed condition table with stable simulation IDs and colors.

    Pass conditions extracted from the source generator for a particular study.
    The default preserves the existing plotting API and archived scan ordering.
    """
    conditions = simulation_conditions() if conditions is None else conditions
    if not conditions:
        raise ValueError("Simulation conditions must not be empty")
    num_flows = max(condition.boundary_flow_key for condition in conditions) + 1
    records = [
        {
            "sim_key": condition.sim_id,
            "pattern_key": condition.pattern_key,
            "pattern_name": pattern_idx2name[condition.pattern_key],
            "boundary_flow_key": condition.boundary_flow_key,
            "boundary_flow": condition.boundary_flow,
            "color": simulation_color(
                condition.pattern_key, condition.boundary_flow_key + 1, num_flows
            ),
        }
        for condition in conditions
    ]
    return (
        pd.DataFrame.from_records(records)
        .set_index("sim_key")
        .astype(
            {
                "pattern_key": "int64",
                "boundary_flow_key": "int64",
                "boundary_flow": "float64",
                "color": str,
            }
        )
    )


def plot_boundary_flux(df: pd.DataFrame, results_dir: Path) -> None:
    """Save distinct absolute flow magnitudes without opening a graphics window."""
    from matplotlib import pyplot as plt

    results_dir.mkdir(parents=True, exist_ok=True)
    flows = df.drop_duplicates("boundary_flow_key").sort_values("boundary_flow_key")
    figure, axes = plt.subplots(figsize=(5, 5), dpi=300, layout="constrained")
    try:
        axes.plot(
            np.abs(flows.boundary_flow.to_numpy(dtype=float)),
            marker="o",
            markersize=10,
            markeredgecolor="black",
            label="simulation",
        )
        axes.set_ylim(bottom=0)
        axes.set_xlabel("simulation index", fontsize=15, fontweight="bold")
        axes.set_ylabel("absolute boundary flow", fontsize=15, fontweight="bold")
        axes.legend()
        figure.savefig(results_dir / "boundary_flux.png", bbox_inches="tight")
    finally:
        plt.close(figure)


def plot_colors(df: pd.DataFrame, results_dir: Path) -> None:
    """Save the pattern/flow color mapping and close the figure reliably."""
    from matplotlib import pyplot as plt

    results_dir.mkdir(parents=True, exist_ok=True)
    figure, axes = plt.subplots(figsize=(5, 5), dpi=300, layout="constrained")
    try:
        axes.scatter(
            df.pattern_key,
            df.boundary_flow_key,
            c=df.color.tolist(),
            edgecolors="black",
            s=100,
        )
        axes.set_ylim(bottom=0)
        axes.set_xlabel("zonation pattern", fontsize=15, fontweight="bold")
        axes.set_ylabel("boundary flow", fontsize=15, fontweight="bold")
        figure.savefig(results_dir / "boundary_color.png", bbox_inches="tight")
    finally:
        plt.close(figure)


def main(argv: list[str] | None = None) -> int:
    """Inspect both studies without writing to the source or loading HDF arrays."""
    from research.spt.cache import cached_catalog, validate_roots
    from research.spt.catalog import CatalogValue, discover_runs, inspect_run
    from research.spt.cli import selected_studies, selection_parser, write_catalog
    from research.spt.conditions import conditions_from_source

    parser = selection_parser(__doc__ or "Catalogue SPT data")
    parser.add_argument(
        "--plots",
        action="store_true",
        help="Save condition plots without showing windows",
    )
    parser.add_argument(
        "--excel", action="store_true", help="Also export the catalogue as Excel"
    )
    parser.add_argument(
        "--dry-run", action="store_true", help="Print metadata without writing reports"
    )
    parser.add_argument(
        "--refresh",
        action="store_true",
        help="Reread descriptors instead of reusing current catalogue metadata",
    )
    args = parser.parse_args(argv)
    rows: list[dict[str, CatalogValue]] = []
    try:
        source, target = validate_roots(args.input_root, args.output_root)
        runs = discover_runs(source, selected_studies(args.study), tuple(args.runs))
        for run in runs:
            try:
                row = None if args.refresh else cached_catalog(run, target)
                rows.append(
                    {
                        **(inspect_run(run).row if row is None else row),
                        "status": "available",
                    }
                )
            except (OSError, ValueError, KeyError) as error:
                rows.append(
                    {
                        "study": run.study,
                        "run_id": run.run_id,
                        "status": "failed",
                        "error": str(error),
                    }
                )
        if args.study in {"both", "simulations"} and not args.runs:
            conditions, _ = conditions_from_source(source)
            present = {run.run_id for run in runs if run.study == "simulations_sixth"}
            for condition in conditions:
                if condition.sim_id not in present:
                    rows.append(
                        {
                            "study": "simulations_sixth",
                            "run_id": condition.sim_id,
                            "status": "missing",
                        }
                    )
        data = pd.DataFrame(rows)
        print(
            data.to_string(
                index=False,
                columns=[
                    name
                    for name in (
                        "study",
                        "run_id",
                        "num_cells",
                        "num_steps",
                        "tend_s",
                        "complete",
                        "status",
                    )
                    if name in data
                ],
            )
        )
        if not args.dry_run:
            write_catalog(data, target, excel=args.excel)
            if args.plots and args.study in {"both", "simulations"}:
                conditions, _ = conditions_from_source(source)
                frame = simulation_conditions_df(conditions)
                plot_boundary_flux(frame, target / "figures")
                plot_colors(frame, target / "figures")
    except (OSError, ValueError) as error:
        print(f"SPT information: {error}", file=sys.stderr)
        return 1
    return int(any(row.get("status") in {"failed", "missing"} for row in rows))


if __name__ == "__main__":
    raise SystemExit(main())
