"""Name every directory and record of the SPT results tree in one place.

The results tree separates deliverables from intermediate data::

    results/
    ├── index.html      offline report; open this first
    ├── README.md       this layout, written with every report
    ├── manuscript/     manuscript figures, tissue panels and stable exports
    ├── analysis/       quantitative study analysis: tables, NetCDF and figures
    ├── figures/        per-run geometry time courses and study overviews
    ├── meshes/         reconstructed lobules and mixed-sector histories
    ├── catalog/        run catalogue and simulation-condition figures
    ├── logs/           execution reports and per-stage catalogues
    └── cache/          intermediate data, rebuilt by rerunning the pipeline
"""

from pathlib import Path

CATALOG = "catalog"
LOGS = "logs"
MESHES = "meshes"
MANUSCRIPT = "manuscript"
EXPORTS = "exports"
REPORT = "index.html"
README = "README.md"

# Per-run or per-study artifact categories, relative to the results root.
CATEGORIES: dict[str, tuple[str, ...]] = {
    "cells": ("cache", "cells"),
    "analysis-data": ("cache", "analysis"),
    "plot-data": ("cache", "plots"),
    "plot-depth": ("cache", "depth"),
    "analysis": ("analysis",),
    "figures": ("figures", "runs"),
    "overviews": ("figures", "overviews"),
    "meshes": (MESHES,),
}
PREVIEWS = ("cache", "previews")


def category_root(root: Path, category: str) -> Path:
    """Return the directory holding one artifact category for every study."""
    return root.joinpath(*CATEGORIES[category])


def artifact_directory(root: Path, category: str, study: str, run_id: str) -> Path:
    """Return the stable manifest directory of one run within a category."""
    return category_root(root, category) / study / run_id


def catalog_path(root: Path, suffix: str = ".tsv") -> Path:
    """Return the run catalogue written by the catalogue stage."""
    return root / CATALOG / f"runs{suffix}"


def log_path(root: Path, name: str) -> Path:
    """Return one execution record below the logs directory."""
    return root / LOGS / name


def manuscript_directory(root: Path, stem: str) -> Path:
    """Return the manifest directory of one manuscript figure."""
    return root / MANUSCRIPT / stem


def exports_directory(root: Path) -> Path:
    """Return the directory of stable manuscript filenames."""
    return root / MANUSCRIPT / EXPORTS


README_TEXT = """# SPT results

Open `index.html` for the offline report of every figure, table and record.
Keep this directory together; the report uses relative links.

| Directory | Contents |
| --- | --- |
| `manuscript/exports/` | Stable manuscript filenames: `figure-04.pdf` to `figure-07.pdf`, `figure-s1-convergence.pdf`, compact variants and the linear-increase tissue panels, each with caption and provenance |
| `manuscript/<figure>/` | Versioned generation of each manuscript figure |
| `analysis/<study>/` | Quantitative analysis of homogeneous lobules: `analysis.nc`, time-course and spatial-profile tables, figures and provenance |
| `figures/runs/<study>/<run>/` | Geometry time course of one sixth, lobule or mixed-sector history |
| `figures/overviews/<study>/` | Final-state overview of all runs of a study and geometry |
| `meshes/<study>/<run>/` | Reconstructed lobules and mixed-sector XDMF/HDF histories with comparison previews |
| `catalog/` | `runs.tsv` run catalogue (`runs.xlsx` with `--excel`) and simulation-condition figures |
| `logs/` | Pipeline, processing, analysis, plotting and manuscript reports with per-stage catalogues |
| `cache/` | Intermediate numerical caches, plot samples and report previews |

Studies: `simulations_*` are the 48 zonation-pattern and dose simulations,
`convergence_*` the mesh-convergence runs. `*_sixth` is the simulated sixth of a
lobule, `*_lobulus` the reconstructed homogeneous lobule and `*_combinations`
six different sixths assembled for visual comparison.

Each artifact directory holds a `manifest.json` that names its current
`generation-<id>/` folder. Only the current generation is kept; the empty
`.lock` file lets concurrent runs coordinate. Rerun the pipeline to update or
repair results: unchanged results are reused and leftovers of interrupted runs
are removed.
"""
