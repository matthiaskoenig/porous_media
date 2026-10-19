"""Shared command-line selection and tabular catalogue exports."""

import argparse
from pathlib import Path
from typing import cast

import pandas as pd

from research import RESULTS_DIR
from research.spt import data_spt_dir
from research.spt.artifacts import atomic_output
from research.spt.catalog import STUDIES, Study


def selection_parser(description: str) -> argparse.ArgumentParser:
    """Provide explicit input/output roots and stable study/run selections."""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--input-root", type=Path, default=data_spt_dir)
    parser.add_argument("--output-root", type=Path, default=RESULTS_DIR / "spt")
    parser.add_argument(
        "--study", choices=("both", "simulations", "convergence"), default="both"
    )
    parser.add_argument("--runs", nargs="+", default=[])
    return parser


def selected_studies(selection: str) -> tuple[Study, ...]:
    """Resolve user-facing study names to the supplied directory names."""
    if selection == "both":
        return STUDIES
    if selection not in {"simulations", "convergence"}:
        raise ValueError(f"Unknown study: {selection}")
    return (cast(Study, f"{selection}_sixth"),)


def write_catalog(
    data: pd.DataFrame, output_root: Path, *, excel: bool = False, name: str = "catalog"
) -> None:
    """Write selected-run metadata locally, with Excel as an optional export."""
    for suffix in (".tsv", ".xlsx") if excel else (".tsv",):
        with atomic_output(output_root / f"{name}{suffix}") as temporary:
            if suffix == ".tsv":
                data.to_csv(temporary, sep="\t", index=False)
            else:
                data.to_excel(temporary, index=False)
