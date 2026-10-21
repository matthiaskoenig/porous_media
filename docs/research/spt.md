# SPT simulations and convergence

The repository-only SPT workflow catalogues supplied XDMF/HDF results,
reconstructs lobules, serializes quantitative analysis, and generates numerical
and geometry figures. It requires access to study data, but does not need FEBio,
VTK rendering, or a graphics session. The reusable library remains independent
of these study workflows. This analyzes existing solver results; it does not run
new simulations.

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
uv run --no-sync python -m research.spt.data_processing --step all --output-root results/spt
```

To run the complete project from catalogue through analysis and figures, use:

```bash
uv run --no-sync python -m research.spt.pipeline --output-root results/spt
```

The [complete workflow](#run-and-verify-the-complete-project) below describes
resuming, inspecting stages, and the resulting project index.

On PowerShell, set `$env:POROUS_MEDIA_SPT_DATA` and `$env:MPLBACKEND` instead.
Both commands also accept `--input-root /data/koeniglab/spt` explicitly. The input
default is `POROUS_MEDIA_SPT_DATA`, falling back to the checkout's `data/spt`.
The output default is `RESULTS_DIR / "spt"`, normally `results/spt`, on local
storage. Overlapping input/output trees and symlinked output subdirectories are
rejected.
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

`data_processing.py` uses existing XDMF/HDF results directly. Its default cache
step reads the source sixths. Explicit reconstruction steps prepare full lobules
for subsequent analysis, as described below. VTK conversion remains separate.

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

`plots_simulations.py` and `plots_convergence.py` are entry points to the shared
[quantitative analysis](#time-courses-spatial-profiles-and-convergence) below.
They use actual run IDs and local homogeneous-lobule caches. Geometry plotting
uses its own shared renderer. Numerical cache generation does not invoke either
plotting workflow or import PyVista.

Lobule reconstruction is an explicit processing step described below. VTK
conversion remains a separate helper for studies without XDMF results.
PyVista snapshots and video export require a working rendering setup; the sector
comparison figures generated below use Matplotlib and run headlessly.
For a portable numerical example, see the [library quickstart](../quickstart.md).

## Reconstruct lobules and sector comparisons

Create the modified geometries after inspecting the supplied sixths. This step
reconstructs existing simulation results by symmetry; it does not run FEBio.
The complete workflow then builds numerical caches from the homogeneous lobules:

```bash
uv run --no-sync python -m research.spt.data_processing --step all --workers 2
```

Individual steps can be run or resumed independently:

```bash
# Reconstruct every simulation and convergence mesh, preserving native times.
uv run --no-sync python -m research.spt.data_processing --step lobulus

# Compose the mixed convergence view and eight comparisons of zonation patterns.
# This also saves PNG and PDF figures alongside the time series.
uv run --no-sync python -m research.spt.data_processing --step combinations

# Build both sets of geometries and figures, without numerical cache creation.
uv run --no-sync python -m research.spt.data_processing --step geometry

# Analyze reconstructed lobules using the same scalar profiles and cache policy.
uv run --no-sync python -m research.spt.data_processing \
  --step cache --geometry lobulus

# Inspect reconstruction status without writing or starting workers.
uv run --no-sync python -m research.spt.data_processing --step geometry --status
```

The existing default, `--step cache --geometry sixth`, remains available for
inspecting source-sector results. Subsequent full-lobule analysis should use
`--geometry lobulus` or `--step all`. Mixed comparisons are excluded from this
quantitative analysis path. A selected reconstruction can use `--step lobulus
--runs sim001`; combinations require complete studies and reject `--runs`.
`--study simulations` and `--study convergence` work for every step.

The reconstruction writes:

| Study | Contents | Purpose |
| --- | --- | --- |
| `simulations_lobulus` | All 48 simulation lobules | Quantitative analysis and visualization |
| `convergence_lobulus` | All seven convergence lobules | Quantitative analysis and visualization |
| `simulations_combinations` | `dose_01` to `dose_08`, one sector per pattern | Visualization |
| `convergence_combinations` | `convergence_mixed`, one sector per resolution | Visualization |

Each lobule places six copies of its source sector counterclockwise around the
origin, reflecting alternating sectors to pair matching radial boundaries.
Connectivity permutations preserve each source cell's orientation. Scalar fields
are copied, Cartesian vectors are transformed with Q, and full stress/strain
tensors with Q T Q.T. The Cartesian scalar fields `position_x`, `position_y`, and
`position_z` transform together. Unsupported component layouts fail explicitly.
Cell volume is preserved per cell; a homogeneous lobule has six times its
source's volume and the same volume-weighted means and standard deviations.

Seam points remain separate. This preserves different values on adjoining
patterns and permits different convergence resolutions without inventing node
correspondence or averaging discontinuities. `sector_id` identifies sectors
0 through 5, and `source_cell_id` preserves the original cell identity within
each source. These are reconstructed visualization/analysis geometries, not
conforming meshes intended as new solver inputs.

The mixed convergence view uses `06250`, `09375`, `12500`, `18750`, `25000`, and
`37250` in that order. `50000` still receives its own homogeneous lobule. The
finest input ends at 36,590 s, so the mixed history ends there. Native output times
are offset between resolutions. The mixed history uses the first sector's
native time grid inside the common interval, including the common endpoints.
Continuous scalar/vector/tensor fields interpolate linearly. Integer fields,
`cell_type`, `rr_zonation_pattern`, `rr_necrosis`, and `gauss_points_TPM` retain
their previous sampled value. No field is extrapolated. Homogeneous lobules
retain every original time sample and do not interpolate any field.

Each dose comparison combines patterns 0 through 5 from the same boundary-flow
condition. For example, dose 1 uses `sim001`, `sim009`, `sim017`, `sim025`,
`sim033`, and `sim041`. Dose numbers identify the generator's boundary-flow
ordering; flow units remain unspecified. Each comparison includes the complete
common history. Its figure shows initial protein zonation and final necrosis,
with numbered sectors and a pattern legend. The convergence figure shows the
six meshes and final toxicant values. Figures render headlessly with Matplotlib,
without requiring a VTK display backend.

Geometry generations are stored separately from numerical cache generations:

```text
results/spt/<study>/xdmf/<run_id>/
├── manifest.json
└── generation-<identifier>/
    ├── <run_id>.xdmf
    ├── <run_id>.h5
    ├── <run_id>_limits.json
    ├── provenance.json
    ├── model.feb          # Homogeneous lobules, when a source model exists.
    ├── <run_id>.png       # Mixed comparisons only.
    └── <run_id>.pdf       # Mixed comparisons only.
```

XDMF references its colocated HDF sidecar by relative filename. Copying the pair
is sufficient to open it elsewhere; it does not depend on network HDF files.
Each worker stages its input descriptors and HDF sidecars on local storage with
sequential copies, avoiding repeated small reads over the network filesystem.
Temporary copies are removed after readers close, including on failure. Allow
space for one source history per lobule worker or six per mixed-view worker,
in addition to generated outputs. Original inputs remain the provenance source
and are checked for changes before publication.
The HDF writer uses gzip compression and shares identical arrays within the
sidecar to reduce repeated scalar and time-invariant storage. Source fingerprints,
conditions, field-transform policy, and temporal alignment are recorded in the
manifest and provenance file. Only complete generations become visible through
the manifest. Failed reconstruction or figure generation preserves the previous
generation, and unchanged runs skip array reads and worker startup.
`geometry-report.json` describes the latest selected reconstruction invocation.
`geometry-catalog.tsv` lists successful outputs, including descriptor and figure
paths, so readers do not need to navigate generation directories manually.

Resolve published paths rather than globbing temporary or retired generations:

```python
from pathlib import Path

from research.spt.cache import load_summary
from research.spt.geometry import geometry_path, lobulus_sources

raw = Path("/data/koeniglab/spt")
output = Path("results/spt")
source = lobulus_sources(raw, output, run_ids=("sim001",))[0]
summary = load_summary(source, output)
lobulus_xdmf = source.xdmf
comparison_xdmf = geometry_path(
    output / "simulations_combinations" / "xdmf" / "dose_01"
)
```

Numerical lobulus caches use the same `runs/<run_id>/manifest.json` layout as
source-sector caches, under their respective `*_lobulus` study directories.
Use published lobules for quantitative analysis and sector comparisons for
visualization. Mixed-sector histories never enter whole-lobule statistics.

## Plot sixths, lobules and combined geometries

After reconstruction, create all geometry figures with:

```bash
uv run --no-sync python -m research.spt.plot_geometry --workers 2
```

The input, data and output roots default to the same locations used above.
`--data-root results/spt` selects the published reconstructed histories;
`--output-root results/spt` selects local plot storage. The raw share is read
only. Plotting does not change reconstructed histories or numerical caches.
These repository commands also work independently:

```bash
uv run --no-sync python -m research.spt.plots_simulations_geometry
uv run --no-sync python -m research.spt.plots_convergence_geometry

# Render a selected source sixth, its lobule and the same-dose comparison.
uv run --no-sync python -m research.spt.plot_geometry \
  --study simulations --runs sim048 dose_08 --geometry all

# Inspect pending work without writing files or reading numerical arrays.
uv run --no-sync python -m research.spt.plot_geometry --status
```

`--geometry` accepts `all`, `sixth`, `lobulus`, or `combined`. The complete current
dataset produces 119 time-course figures: 55 source sixths, 55 homogeneous
lobules, eight dose comparisons and one mixed convergence history. Each of the
six study/geometry groups also gets a final-state overview and an overview of
distinct meshes. Identical simulation meshes appear once in the mesh overview.
Simulation overviews follow the generator's pattern and dose ordering.
Convergence histories additionally export a separate `mesh.svg`, `mesh.pdf` and
`mesh.png`. These neutral mesh-only views show cell boundaries independently of
the scalar fields and appear in the report's **Meshes only** section.

Default figures preserve the manuscript's five snapshots: indices 1, 3, 5, 7
and 9 of ten equally spaced times in the common interval. Explicit `--times`
and other snapshot counts override this selection. The five columns retain
extracellular substrate and product, protein, toxic compound and necrosis.
The canonical layers in `research/spt/fields.py` preserve `magma` for
concentrations, `hot` for protein and `binary` for necrosis, with the original
labels and horizontal scalar bars in every panel. Simulation cell-edge strokes
use 18% opacity while scalar colors remain fully opaque. Convergence scalar
panels omit internal mesh lines and show only the domain outline, including
holes. Homogeneous lobule outlines omit artificial sector seams; independent
combined sectors retain their own outlines. Mesh-only figures use a white
background and 45% opacity cell-edge strokes.
Each panel occupies 4 by 4 inches, matching the original 600 pixel panels at
150 DPI. Scalar bars retain five labels, bold text and one decimal place.
Final-state overviews use the same panel sizes and per-panel scalar bars.
Small ranges use scientific notation when fixed decimals would repeat labels,
so nonzero product concentrations do not appear as five identical zeros.

Times share a common interval within each selected study family. Simulation
plots therefore retain 24 h even when convergence ends at 36,590 s. Shared color
ranges use validated local numerical-cache limits, then full-history
`_limits.json` metadata beside source descriptors, without scanning network
arrays again. Missing or malformed legacy metadata falls back to sampled
ranges when no validated cache is available. Explicit layer limits take precedence. Combined
views identify each sector. These settings are shared with the historical
`DataLayer` adapter so the two workflows cannot silently choose different colors.

Continuous scalars are linearly interpolated between neighboring native frames.
Necrosis, zonation, cell type and other discrete fields retain their previous
native value. No extrapolation or averaging between cells is performed.
Only requested scalar arrays and neighboring frames are read from HDF; point,
vector, tensor and unused scalar arrays are omitted. A compact local sample cache
is shared by time courses, overviews and optional videos. Plotting never creates
interpolated XDMF/HDF histories.

Choose fields, times and outputs explicitly:

```bash
# Final-state figures for only necrosis and toxicant.
uv run --no-sync python -m research.spt.plot_geometry \
  --fields rr_necrosis 'rr_(T)' --snapshots 1 --formats svg pdf

# Explicit times in seconds must lie within every selected history's interval.
uv run --no-sync python -m research.spt.plot_geometry \
  --study simulations --times 0 21600 43200 86400

# Preserve outliers while showing smaller concentrations on a symlog color scale.
uv run --no-sync python -m research.spt.plot_geometry \
  --study convergence --scale symlog --output-root results/spt/convergence-symlog

# Generate only the two study overviews per geometry.
uv run --no-sync python -m research.spt.plot_geometry --only-overviews

# Separate field/time panels are optional, as are videos and GIFs.
uv run --no-sync python -m research.spt.plot_geometry \
  --study simulations --runs dose_08 --geometry combined \
  --panels --video-frames 200 --fps 12 --gif
```

SVG and PDF preserve vector cell geometry and text. PNG previews default to
150 DPI; `--dpi` also controls optional video resolution. FFmpeg is required only
when videos are requested. Video artists are reused and frames stream directly
to FFmpeg, without writing intermediate PNG panels. GIF conversion is opt-in.
`--no-overviews` omits study overviews. Bounded spawned workers own their files
and rendering state; `--workers 1` runs serially.

Dense SVG time courses can be much larger than compressed vector PDFs. Use
`--formats pdf png` for smaller batch outputs and reserve SVG for figures that
need cell-level vector editing. Selecting fewer fields or snapshots also reduces
rendering time and storage without changing the underlying histories.

The shared Matplotlib renderer runs headlessly without importing PyVista or
needing a graphics session. It supports planar XY triangles/quads and vertical
wedge/hexahedron extrusions. Other geometries fail clearly instead of producing
an incorrect projection; use the library's 3D visualization helpers for those.
Historical `visualize_spt_2d` calls remain available, including generator inputs,
custom layer labels, unit conversions and fixed color limits. Separate panels
and videos now require explicit selection in that helper too.

Outputs and preview caches use the following layout:

```text
results/spt/
├── plot-data/<study>/<run_id>/manifest.json
├── figures/<study>/<run_id>/manifest.json
├── overviews/<study>/manifest.json
├── report-previews/<fingerprint>.webp
├── figures/index-both.html
├── plot-catalog-both.tsv
└── plot-report-both.json
```

Each manifest points to a complete generation. A failed replacement preserves
the previous published figures. Source and output fingerprints enable unchanged
runs to skip XML parsing, numerical reads, rendering and worker startup. Changing
only a rendering option reuses the existing sampled arrays. Missing artifacts
are regenerated; `--force` explicitly resamples and rerenders. Fingerprints use
size and modification time, so changes preserving both require `--force`.
The HTML gallery links previews, vector figures and any requested videos; the
catalogue gives their current paths. The report records errors, sampled color
ranges, elapsed time and the number of selected HDF field arrays read. Files
are named by the selected study (`both`, `simulations`, or `convergence`) and
describe the latest invocation, including an intentionally partial selection.

`--scale symlog` is an explicit alternative for concentration ranges containing
large outliers. It preserves every value, with a linear region within one
millionth of the maximum absolute display limit and logarithmic spacing outside
that region. Figures label the transformed scale; necrosis stays linear.
The current `convergence_50000` has sampled substrate values above 2,000 mM,
whereas the six finer histories stay near 1 mM or below. Its product values
likewise stretch the shared linear range. Keep those source values visible;
use a separate symlog output root when inspecting finer spatial differences.

## Time courses, spatial profiles and convergence

After numerical cache generation, run both study analyses or select one:

```bash
uv run --no-sync python -m research.spt.analysis
uv run --no-sync python -m research.spt.plots_simulations
uv run --no-sync python -m research.spt.plots_convergence

# Change spatial resolution or explicit comparison times, in seconds.
uv run --no-sync python -m research.spt.analysis \
  --bins 30 --spatial-times 0 21600 36000 --time-points 101

# Analyze an explicit subset, retaining its source-derived condition labels.
uv run --no-sync python -m research.spt.analysis \
  --study simulations --runs sim001 sim009 --workers 1
```

The default `--geometry lobulus` analyzes all 48 simulation lobules and seven
convergence lobules. `--geometry sixth` explicitly analyzes source-sector
numerical caches instead. Mixed geometries are excluded. `--data-root` selects
the numerical-cache root; `--output-root` selects local analysis storage.
Neither the raw share nor reconstructed simulation histories are modified.

Manuscript time-course figures use matched-time cell-volume-weighted means.
Simulation panels retain the six pattern palettes (`Greys`, `Blues`, `Purples`,
`Reds`, `Oranges`, `Greens`), dose shades, black-edged circular markers and
population SD error bars. `--no-spread` explicitly omits the bars; `--spread`
also enables them for convergence. SD describes cell variation, not uncertainty
in a mean or a confidence interval. Convergence retains a coarse-to-fine
`viridis` palette, plain lines, the 1.5 h reference line and the legend above
the panels. Actual mesh IDs replace obsolete IDs from the archived scripts.
Colors remain stable when selecting a subset of doses or meshes.

Panels retain bold 12 point labels, 2.7 by 2.5 inch dimensions, the original
axis labels and PP/PV ticks, and 300 DPI PNG exports. The native histories,
including every original output and each run's actual endpoint, remain available
in `native-timecourses.tsv` and the per-run `native.nc` caches.

For comparisons, 51 equally spaced times cover the common interval separately
within each selected study family. Continuous **cell values and cell volumes
are interpolated before weighted reduction**. This differs from interpolating
summary means when volumes change. Discrete fields retain their previous native
value, and no time is extrapolated. Source necrosis fractions are preserved;
figures display their volume-weighted percentage while serialized values remain
fractions. The cell data are read selectively from local NetCDF caches, without
loading complete histories, point arrays, or network HDF data.

Spatial plots default to the first common time, the historical reference time
(10 h for simulations and 1.5 h for convergence), and the last common time.
If the reference lies outside the common interval, its midpoint is used instead.
Figures show individual cell scatters rather than binned curves, including
valid fractional necrosis. For quantitative comparisons, twenty equal bins cover the source's normalized `rr_position`, from
periportal 0 to perivenous 1. Each bin stores volume-weighted mean, population
SD, minimum, maximum, cell count, and total volume. Position 1 belongs to the
last bin. Empty or zero-volume bins retain missing statistics and zero volume;
they are never filled by spatial interpolation. Raw cell snapshots are also
serialized so later analysis can choose a different reduction.

The initial-properties figure always uses each run's **first native frame**,
even when explicit spatial comparison times start later. Simulation diagnostics
show the highest selected dose per pattern; convergence diagnostics show every
mesh. Diagnostics retain the manuscript's ten-bin position count histogram and
individual cell scatters. Display conversions follow the archived script:
mesh cell volume uses m³ to nL (`1e12`), and fluid/solid volume uses L to nL
(`1e9`). Pressure retains the original `Pressure [?]` label. Serialized arrays
remain in source units; caches do not provide independent unit metadata.
The former use of the last frame for initial diagnostics remains corrected.

Simulation figures compare doses by pattern and include matched-time dose
responses. Convergence figures use the actual mesh IDs and lobulus cell counts.
Resolution panels retain a linear element-count axis, black connecting line
and mesh-colored markers at 1.5 h. If that time is absent from explicitly
selected snapshots, the middle selected snapshot is used and the report names
its time. The mesh with the largest cell count supplies the reference. Exported tables
report absolute and relative differences of whole-lobule means at matched times.
A zero reference makes relative differences undefined, stored as missing values.
Spatial comparisons report RMS differences of bin means, weighted by reference
bin volume, with the fraction of reference volume covered by mutually populated
bins. These are differences from the finest **available** mesh, not estimates
of error against a continuum solution or an inferred convergence order.

Convergence includes supplemental symlog figures to expose smaller values beside
the current coarsest mesh's large substrate/product outliers. `--linear-only`
omits those supplements. All values remain in the serialized results. Default
exports are SVG/PDF/PNG; `--formats`, `--dpi`, and `--fields` limit figure output.

### Reusable analysis results

Each run publishes a compact numerical generation under
`analysis-data/<study>/<run_id>/`:

| File | Contents |
| --- | --- |
| `native.nc` | Native-time mean and population SD for selected scalars |
| `aligned.nc` | Matched-time statistics, reduced after cell interpolation |
| `profiles.nc` | Spatial statistics, counts, bin boundaries, and volumes |
| `cells.nc` | Selected cell fields at requested spatial times |
| `initial.nc` | Binned diagnostics from the first native frame |

Each study publishes an independent generation under `analysis/<study>/`,
containing `analysis.nc`, native/matched time-course TSV tables, spatial-profile
TSV tables, figures, a local HTML gallery, and provenance. Convergence also
includes `convergence-means.tsv` and `convergence-profiles.tsv`. The study NetCDF
combines matched statistics and profiles with run, dose, pattern, cell-count,
and native-end-time coordinates; native histories with differing time axes are
preserved in the native TSV table and per-run files. No pickle files are used.
The study tables and NetCDF files can be moved and read without the raw share.

`analysis/index.html` links the selected study galleries.
`analysis-catalog-<study>.tsv` supplies current result paths, and
`analysis-report-<study>.json` records errors, common windows, elapsed time, and
the number of native frames read. These describe the latest selection. Unchanged
invocations inspect fingerprints without opening numerical arrays or starting
workers. Plot styling or format changes reuse numerical reductions. Missing
exports rebuild, failed publication retains previous figures, and `--force`
explicitly regenerates the selected analysis. `--status` writes nothing.

Load a validated per-run result from a checkout:

```python
from pathlib import Path

from research.spt.analysis_data import load_analysis
from research.spt.plot_data import plot_sources

raw = Path("/data/koeniglab/spt")
output = Path("results/spt")
source = plot_sources(
    raw, output, study="simulations", geometry="lobulus", run_ids=("sim001",)
)[0]
with load_analysis(source, output, "profiles") as profiles:
    print(profiles.sel(field="rr_necrosis"))
```

### Run and verify the complete project

```bash
uv run --no-sync python -m research.spt.pipeline --status
uv run --no-sync python -m research.spt.pipeline --workers 2
```

The pipeline runs these stages in order:

1. Catalogue actual inputs and plot simulation conditions.
2. Reconstruct homogeneous lobules and mixed visualization geometries.
3. Build or resume the homogeneous-lobule numerical caches.
4. Serialize quantitative analysis and create its figures.
5. Create or resume sixth, lobulus, and combined geometry figures.

It checks every selected result, stops dependent work after a failure, and writes
`pipeline-report.json` plus the final HTML report at the output root. The JSON names
each stage and records its outcome and elapsed time. Repeating the command is a
supported resume operation. Existing numerical and geometry outputs are reused
only when their provenance and output fingerprints match; they are not rebuilt
merely to verify the pipeline.
Short source histories are listed separately in `truncated_runs`; a successful
pipeline does not imply that every source solver run reached its model duration.

`--study` selects a study, and `--runs` selects an explicit subset. A subset
omits mixed comparisons, which require a complete study. `--no-geometry-figures`
skips the last stage when only numerical results are needed. `--formats pdf png`
reduces export size, while `--geometry-snapshots` changes the geometry figure
panel count. Inspecting a fresh output root reports downstream stages as pending
without creating directories or loading HDF arrays.

Portable regressions run every stage from fresh synthetic raw files, including
mixed convergence meshes, installed-library imports, serialized results and
figures. Run them with `uv run --no-sync pytest -q research/tests`. Validation of
the private current dataset uses the same pipeline against the read-only share.

### Shared workflow code

SPT modules share field definitions in `fields.py`, filesystem publication and
locking in `artifacts.py`, and bounded job execution in `execution.py`. Numerical
caches, reconstructed geometries, analysis, and plotting use those same helpers.
Serial and parallel batches isolate job failures consistently, and empty batches
or a single pending history do not start worker processes. Field ordering,
manifest schemas, existing cache generations, and command names remain compatible.

The four `plots_*` command modules delegate to the shared analysis or geometry
workflow. Historical plotting adapters still serve other research projects and
delegate to the current renderer. These helpers remain repository-only research
code, outside the public library and release artifacts.

### Publication figures

Create the updated Figures 4-7 of the
[SPT manuscript](https://doi.org/10.1101/2024.03.26.586870) and a convergence
supplement after the complete pipeline and geometry plotting have finished:

```bash
MPLBACKEND=Agg PYVISTA_OFF_SCREEN=true uv run --no-sync \
  python -m research.spt.publication --data-root results/spt
```

This step reads validated local NetCDF summaries, selected cell frames and cached
projected meshes. It does not require the raw network drive, rerun reductions or
alter existing plots. Unchanged figures use metadata-only cache checks. Missing
exports are repaired; `--force` regenerates selected publication figures.

| Figure | Data and layout | Restored annotations |
| --- | --- | --- |
| 4 | `sim014`, linear increase, D6; five fields by five times | Pattern/dose heading, field headers, time labels and arrow; one bottom colorbar per field |
| 5 | `sim006`, `sim014`, `sim022`, `sim030`, `sim038`, `sim046` at D6; two columns by three pattern blocks | Consistent A-F pattern labels, field headers, time labels and arrows; one shared set of five bottom colorbars |
| 6 | All 48 homogeneous lobules; six pattern rows and five fields | Consistent A-F pattern labels, field headers and D1-D8 legend |
| 7 | Necrosis at 24 h; six pattern rows by eight doses | Consistent A-F pattern labels, D1-D8 columns, inflow arrow and one large vertical necrosis colorbar on the right |
| S1 | Six source-sector meshes, five common-time mixed-sixth readouts and native lobule time courses | Mesh IDs and cell counts, matched boundary/curve colors, A-D panels, exact map time, linear/symlog labels, common-window boundary and truncated-history markers |

Figures 4 and 5 explicitly sample **0, 6, 12, 18 and 24 h**, independently of the
historical overview snapshot selection. Labels are derived from the actual times;
a nonzero overview frame is never relabeled as zero. Continuous cell fields are
interpolated locally and necrosis holds its preceding native frame. Shared scalar
limits use full histories of the entire simulation cohort, preserving the mesh
palettes. Cell-edge widths scale with panel size while retaining their existing
opacity. Figure 6 reuses the existing volume-weighted means, population SD,
markers and pattern/dose colors. SD represents spatial heterogeneity, rather than
uncertainty across replicate simulations.

Figures 4 and 5 use five shared horizontal colorbars, one per field. In Figure 5,
that single bottom set serves all six pattern blocks. Figure 7 uses one large
vertical necrosis colorbar beside the complete grid, with 0 marked alive and 1
marked death. Individual map panels have no colorbars. The field palettes and
full-cohort limits remain unchanged. Geometry axes use their cropped domain
extents, including white living tissue; smaller canvases and closely packed grids
reduce whitespace around the meshes.

All manuscript figures use dose labels **D1-D8**. The representative dose in
Figures 4 and 5 is **D6**. Pattern names and lettering match across Figures 4-7:
**A Constant**, **B Linear increase**, **C Sharp perivenous**, **D Linear decrease**,
**E Sharp periportal** and **F Random**. Figures 6 and 7 share the same row-label
format and typography; Figure 5 uses those labels as its pattern-block headings.

The convergence supplement retains every available native time point and actual
end time. The `06250` and `09375` histories stop at 10.16 and 21.58 h; no curve is
extended to 24 h. The common comparison window ends at 10.16 h. Both linear and
symlog views use only the six approved meshes `06250`, `09375`, `12500`, `18750`,
`25000` and `37250`, ordered fine to coarse. Necrosis axes stay linear. Each
resolution has a distinct qualitative color, shared by its time courses,
source-mesh boundary and mixed-sector boundaries. Source sixths are selected from
local lobule meshes, so producing this figure does not reopen the original network
meshes. The mixed view uses those same six resolutions.
Panel B compares substrate S, product P, protein, toxic compound T and necrosis
at **5400 s (1.5 h)**, the existing convergence reference time. All six sectors
use exactly this time, native-frame continuous interpolation and preceding-frame
necrosis states. The scalar values come from local homogeneous-lobule cell caches,
with composition order and cell counts checked before assembly. These maps show
sector outlines without internal mesh edges. Each field uses the canonical palette
and full-history limits shared by the six meshes; necrosis uses the fixed range
0-1. Use
`--comparison-time SECONDS` to select another time within every native history;
times outside the common interval are rejected. Dashed time-course lines mark
the map time; dotted lines mark the common interval end.

Stable manuscript filenames are available under `publication/exports/`, for
example `figure-04.pdf` through `figure-07.pdf` and `figure-s1-convergence.pdf`.
These files share storage with validated generations. Regenerating a figure
refreshes its stable exports, caption and provenance from the complete generation.
Outputs live under `publication/<figure>/manifest.json` in complete generations.
Each generation contains PNG, SVG and PDF exports, `caption.md` and
`provenance.json`, including run conditions, snapshot times and native end times.
Dense map panels are rasterized at the chosen export DPI to keep the composite
files practical; annotations and time-course curves stay vector in SVG/PDF.
Each map is rendered as an independent small tile so PDF export never retains
one full-composite raster buffer per panel.
All publication formats are cropped to displayed content with a one-point safety
border. Cropping retains labels, arrows, legends and white living tissue, and
preserves the spacing between panels. Hidden map-axis rectangles do not add
empty margins; embedded map tiles are cropped geometrically, independently of
their scalar colors. SVG/PDF text and time courses remain vector.
The default is 300 DPI. Select `--dpi 600` when higher-resolution maps are needed.
For example, regenerate only Figures 4 and 7 as PDF and PNG:

```bash
uv run --no-sync python -m research.spt.publication \
  --data-root results/spt --figures 4 7 --formats pdf png --dpi 300
```

Create additional, smaller versions of **Figures 4 and 5 only** at **300 DPI** with:

```bash
uv run --no-sync python -m research.spt.publication \
  --data-root results/spt --compact --dpi 300
```

Compact figures use a nominal canvas width of 180 mm, preserving the aspect ratio,
composition, colors and scientific data. Typography, line widths and embedded map
tiles scale with the canvas. PNG resolution remains 300 DPI; SVG/PDF annotations
and curves remain vector. Tight cropping can slightly change the exported width.
Use `--width-mm 200` for another nominal width. Compact exports have independent
caches and stable names such as `figure-04-compact.pdf` and
`figure-05-compact.png`; full-size exports remain available. Figures 6, 7 and S1
use their full-size exports. The publication command retires previously generated
compact versions of those three figures from the export directory and report.
Explicitly requesting them with `--compact` is rejected. Width and variant are
recorded in the provenance. Both sizes of Figures 4 and 5 appear together in the
inspection report, using its small cached previews.

The default representative dose is index 6, matching the original manuscript.
`--dose` changes that selection for Figures 4 and 5 and records the choice in
their provenance and captions. Keep the default when reproducing the manuscript.
The final inspection report includes published outputs in its **Publication
figures** section and uses the same small, cached previews as other figures.

### Final inspection report

Open `results/spt/report.html` or `results/spt/index.html` in a browser after the
pipeline finishes. Both open the same offline report. It gathers all current
published figures: numerical analysis, geometry overviews, individual sixths and
lobules, mixed sectors, supplemental views, and simulation conditions. Exports of
the same figure share one card. Previously generated supplements remain included
when the latest pipeline invocation selects a subset of the data.

Search by run, pattern, dose, or title, and combine section, study family,
geometry, pattern, dose, and scale filters. Click a preview for the zoom viewer.
Left and right arrow keys move through the filtered selection; Escape closes the
viewer. SVG and PDF links open the original exports. Data and provenance links
provide the numerical tables, NetCDF results, generation metadata, and execution
records. Aggregate figures match every pattern and dose they contain.

Refresh only the report after generating additional figures:

```bash
uv run --no-sync python -m research.spt.report --output-root results/spt
```

This reads local manifests and output metadata without accessing the raw share,
opening numerical arrays, or rerendering figures. It lists missing, modified, or
invalid generations as omissions instead of linking incomplete outputs. The
pipeline status refers to the latest pipeline invocation; the report does not
independently validate raw-source freshness. Condition PNGs are included when
present. The gallery uses lazily loaded, lossless WebP previews bounded by
768 by 768 pixels. Opening a figure loads its original resolution; PNG, SVG and
PDF downloads remain unchanged. Previews are cached under `report-previews/`
using relative source paths and file metadata. An unchanged report refresh reads
only metadata and reuses previews without decoding the original images. Missing
or changed previews are regenerated locally without rereading simulation data.

The HTML embeds its styling and interaction code and needs no web server or
internet connection. Images and data retain relative links. Move or share the
whole results directory, including its generation and `report-previews` folders,
to keep those links working. The report generator also recognizes relocated
generation files from their manifests. Generated reports and private results stay
outside Git.
