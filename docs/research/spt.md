# SPT simulations and convergence

The repository-only SPT commands catalogue simulation conditions and prepare local
analysis caches from the supplied XDMF/HDF results. They require access to study
data, but do not need FEBio, VTK rendering, or a graphics session. The reusable
library remains independent of these study workflows.

## Input and output directories

Keep the network dataset separate from local outputs. The current dataset is
available at `/data/koeniglab/spt`, with this structure:

```text
spt/
├── simulations_sixth/
│   ├── material/
│   │   ├── feb_manipulation_v2.py
│   │   └── lobule_zonation_pattern.feb.template
│   └── xdmf/
│       ├── sim001.xdmf
│       ├── sim001.h5
│       ├── sim001_cells.nc
│       ├── sim001_points.nc
│       └── ...
└── convergence_sixth/
    ├── material/
    │   └── convergence_*.feb
    └── xdmf/
        ├── convergence_06250.xdmf
        ├── convergence_06250.h5
        └── ...
```

Use the [development environment](../development.md#environment), or sync just
research dependencies with `uv sync --locked --group research --no-dev`.
Commands run from the checkout root:

```bash
export POROUS_MEDIA_SPT_DATA=/data/koeniglab/spt
export MPLBACKEND=Agg
uv run --no-sync python -m research.spt.information --output-root results/spt
uv run --no-sync python -m research.spt.data_processing --output-root results/spt
```

On PowerShell, set `$env:POROUS_MEDIA_SPT_DATA` and `$env:MPLBACKEND` instead.
Both commands also accept `--input-root /data/koeniglab/spt` explicitly. The input
default is `POROUS_MEDIA_SPT_DATA`, falling back to the checkout's `data/spt`.
The output default is `RESULTS_DIR / "spt"`, normally `results/spt`, on local
storage. Overlapping input/output trees and output symlinks are rejected.
Neither command creates or changes files on the source share.

## Inspect conditions and available runs

`information.py` reads descriptor and model metadata without loading HDF arrays.
After processing, it reuses catalogue metadata when source fingerprints match.
Use `--refresh` to reread descriptors explicitly.
It discovers actual convergence IDs and records cell/point counts, available
scalar fields, time ranges in seconds, and completion relative to the model's
explicit control duration. Models with ambiguous or absent duration metadata
have an unknown completion status.

Simulation IDs follow the pattern/flow ordering in the supplied generator. The
command parses literal assignments without executing the archived Python script.
If no generator is supplied, it labels the repository defaults as unverified.
Flow units remain unspecified when the source does not declare them. Convergence
filename tokens remain identifiers; no mesh-length units are inferred from them.

```bash
uv run --no-sync python -m research.spt.information --dry-run
uv run --no-sync python -m research.spt.information --plots --excel
```

The first command writes nothing. The second writes `catalog.tsv`, optional
`catalog.xlsx`, and two condition figures under `figures/`. Figures are saved
without opening windows and are closed after use. Metadata failures and missing
expected simulations are reported and produce a nonzero exit status.

The dataset inspected on 2026-10-19 contains 48 simulation runs and seven
convergence runs. `convergence_06250` ends at 36,590 s and `convergence_09375` at
77,687.6 s, earlier than the model duration of 86,400 s. These runs remain
available, with `complete=False`. Comparisons must choose a shared time window
or explicitly handle different end times.

## Build or resume local caches

`data_processing.py` uses existing XDMF/HDF results directly. It does not rerun
VTK conversion or construct six copies of each mesh for numerical analysis.

```bash
uv run --no-sync python -m research.spt.data_processing --status
uv run --no-sync python -m research.spt.data_processing --workers 2
uv run --no-sync python -m research.spt.data_processing \
  --study convergence --runs convergence_06250 --workers 1
```

`--status` and `--dry-run` inspect metadata and cache state without writing.
`--study` accepts `both`, `simulations`, or `convergence`. `--runs` takes explicit
run IDs; unknown selections and missing/empty study directories fail clearly.
A complete simulation selection requires every expected run. Use `--runs` for
an intentionally partial study.

The default `analysis` profile retains the five analysis scalars, cell volumes,
position, fluid/solid volumes, pressure, zonation pattern, and cell type. It omits
unused point arrays. `--profile full` retains every scalar cell and point field.
Stored field precision is preserved, with type promotion where needed.

First-time adoption of the supplied NetCDF caches validates time coordinates,
field dimensions, and first/middle/last HDF frames for every selected field.
The manifest labels this as **sample validation**, not a complete comparison of
all values. Choose a stronger policy when needed:

```bash
# Compare every selected cached value against HDF before adoption.
uv run --no-sync python -m research.spt.data_processing --legacy-cache full

# Read selected fields directly from HDF, bypassing supplied NetCDF caches.
uv run --no-sync python -m research.spt.data_processing --legacy-cache ignore

# Rebuild a selected local cache from HDF even if its fingerprints match.
uv run --no-sync python -m research.spt.data_processing --runs sim001 --force
```

Rejected legacy caches fall back to HDF and the report records the reason.
Subsequent runs check input sizes and nanosecond modification times, processing
options, conditions, and output fingerprints. Unchanged runs do not reopen XDMF
or numerical arrays. Changed sources rebuild from HDF, avoiding reuse of stale
legacy caches. Changes that deliberately preserve both file size and modification
time require `--force`; metadata fingerprints are not content hashes.

Two spawned workers are the default. `--workers 1` runs serially. Each worker
owns one run's files and memory. Tune the worker count for the network share
rather than the machine's CPU count. An unchanged batch avoids worker startup.

## Outputs and failure recovery

Each run publishes a complete generation under:

```text
results/spt/<study>/runs/<run_id>/
├── manifest.json
└── generation-<identifier>/
    ├── cells.nc
    ├── points.nc
    ├── summary.nc
    └── limits.json
```

The manifest identifies the active generation and records provenance and the
legacy validation method. Scalar limits are calculated from the selected data
without rereading the source. An OS lock prevents concurrent readers/writers
from replacing a run during a load. Partial writes never replace the active
manifest. A failed run leaves previously completed outputs intact, and the next
invocation resumes other completed runs. Source changes during processing cause
failure and require a retry after the upload completes.

`processing-report.json` records all selected successes and failures, while
`processing-catalog.tsv` records the selected successful metadata. These describe
the latest invocation and do not overwrite the separate information catalogue.
Any failed run gives the processing command a nonzero exit status. Short but
readable histories are reported separately from processing failures.

Cache APIs load data into memory and close file handles:

```python
from pathlib import Path

from research.spt.cache import load_processed_run, load_summary
from research.spt.catalog import discover_runs

source = discover_runs(Path("/data/koeniglab/spt"), run_ids=("sim001",))[0]
output = Path("results/spt")
summary = load_summary(source, output)  # Reads only the compact summary.
cells, points, summary = load_processed_run(source, output)
```

Summaries contain native-time volume-weighted population means and standard
deviations, accumulated in float64. Invalid volumes and nonfinite analysis values
fail explicitly. Per-cell arrays retain the original time axis for spatial
analysis and interpolation. Interpolating a summary is not generally equivalent
to interpolating cell values and changing volumes before reduction.

## Further analysis and visualization

The existing `plots_simulations.py`, `plots_convergence.py`, and geometry plotting
scripts remain study-specific plotting code. Their historical dataset selections,
convergence labels, and output paths need review before using them with new runs.
Use the catalogue's actual IDs and the local cache APIs when updating them.
Numerical cache generation does not invoke these scripts or import PyVista.

Lobule reconstruction and VTK conversion remain explicit research/helper steps
for workflows that need them. Full geometry visualization, snapshots, and video
export can use the original XDMF/HDF pairs and require a working rendering setup.
For a portable numerical example, see the [library quickstart](../quickstart.md).
