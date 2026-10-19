"""Shared command-line selection and atomic research report writing."""

import argparse
import json
import tempfile
from pathlib import Path
from typing import cast

import pandas as pd

from research import RESULTS_DIR
from research.spt import data_spt_dir
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


def atomic_json(path: Path, value: object) -> None:
    """Replace a report only after a complete JSON file has been written."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w", encoding="utf-8", dir=path.parent, prefix=".report-", delete=False
        ) as handle:
            temporary = Path(handle.name)
            json.dump(value, handle, indent=2, allow_nan=False)
            handle.write("\n")
        temporary.replace(path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def write_catalog(
    data: pd.DataFrame, output_root: Path, *, excel: bool = False, name: str = "catalog"
) -> None:
    """Write selected-run metadata locally, with Excel as an optional export."""
    output_root.mkdir(parents=True, exist_ok=True)
    for suffix in (".tsv", ".xlsx") if excel else (".tsv",):
        with tempfile.NamedTemporaryFile(
            dir=output_root, suffix=suffix, delete=False
        ) as handle:
            temporary = Path(handle.name)
        try:
            if suffix == ".tsv":
                data.to_csv(temporary, sep="\t", index=False)
            else:
                data.to_excel(temporary, index=False)
            temporary.replace(output_root / f"{name}{suffix}")
        finally:
            temporary.unlink(missing_ok=True)
