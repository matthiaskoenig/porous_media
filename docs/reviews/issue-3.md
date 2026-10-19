# Issue #3: library and research code review

This review starts from commit `c784dc6` and covers usability, performance,
security, resource handling, packaging, and the boundary between reusable APIs
and project workflows. The review uses bundled data and synthetic fixtures;
private datasets and FEBio installations are not prerequisites.

## Findings and improvements

### High: template evaluation could reach Python internals

`create_feb_files` used a normal Jinja2 environment. A template expression could
access `cycler.__init__.__globals__.os` and call Python functions. The reproducer
called `getcwd()` through that chain. Context strings were also emitted without
XML escaping, and missing parameters silently became empty strings.

The renderer now uses an immutable sandbox, strict undefined values, and XML
escaping. Rendering errors propagate before the affected output is opened.
Tests cover internal attribute access, missing values, and quoted Unicode/XML
content. Templates and supplied Python objects still require review: the
[Jinja sandbox documentation](https://jinja.palletsprojects.com/en/stable/sandbox/)
explains its resource and object-capability limits. This is not an operating
system sandbox or FEBio model validator.

### High: generated filenames could escape the output directory

Simulation identifiers were interpolated directly into paths. With an existing
`outputs/model_x` directory, identifier `x/../../escaped` wrote `escaped.feb`
outside `outputs`. Rendering also replaced existing files without an explicit
overwrite choice and printed entire parameter contexts.

Simulation identifiers and exported image names now require a single portable
filename component. Windows device names and separators are rejected on every
platform. FEBio rendering defaults to exclusive file creation, rejects symlink
outputs, and no longer prints contexts. Image export validates paths before
creating plotters or output directories. Tests check traversal, preserved files,
explicit overwrite, and symlink targets. These checks assume the selected parent
directories remain under the caller's control during writes.

### High: orphaned HDF sidecars could be overwritten

VTK conversion and interpolation checked whether the XDMF output existed, then
opened the associated HDF file in write mode. If only the HDF file existed,
`overwrite=False` still allowed its replacement.

Both converters now protect an existing orphaned sidecar and reject symlink
outputs. Interpolation also rejects input/output paths sharing the conventional
HDF sidecar. Successful regeneration always refreshes the associated limits,
including when an older JSON cache is present. Tests preserve unrelated bytes
and verify refreshed limits after conversion.

### High: zonation produced incorrect positions and ignored cell blocks

The old implementation capped both nearest-boundary distances at one mesh unit.
For boundary centers at 0 and 10, an interior cell at 2 received position 0.5
instead of 0.2. It also used only the first cell-data block while looping over
every topology block, reusing local indices and returning an incomplete field.

Cell centers and labels are now combined across all blocks. Separate spatial
indexes find the true periportal and perivenous distances, and results are split
back into the original block order. Tests cover multiple blocks, length-unit
changes, required boundary labels, valid volumes, and preservation of the input
mesh. Both accepted volume names work when preserving existing fields.

### High: scalar dataset conversion failed on unequal cell blocks

`mesh_datasets_from_xdmf` used `np.dstack` on lists of cell arrays. A synthetic
mesh with two triangle cells and one quadrilateral failed with incompatible
dimensions. Squeezing also removed the item dimension for single-cell inputs.

The converter now concatenates cells in mesh order and fills NumPy arrays
directly, without building a DataFrame and xarray object for each time step.
It preserves `(time, cell)` and `(time, point)` dimensions, promotes changing
numeric types, and rejects changing scalar schemas before writing caches.
Cached and newly generated datasets are fully loaded, with no retained file
handles. Tests cover mixed blocks, one-cell inputs, empty point fields, fractional
values, schema errors, and cache reuse.
Metadata also identifies six-component tensors consistently for points and cells
and records the original array shape.

### Medium: parallel conversion returned large unused results

Directory conversion returned complete datasets from each process even though
the caller only needed cached files. This incurred serialization, interprocess
transfer, and parent memory use. The pool size depended on CPU count without
considering the number of input files.

Workers now write caches and return no datasets. The public function accepts
`max_workers`, bounds its default by available CPUs and file count, handles empty
directories without a pool, and runs serially with one worker. Serial and
two-process tests produce the expected cache pairs. Workers use spawn rather
than forking an application that may already have native library threads.

### Medium: imports changed application state and loaded rendering code

Importing mesh zonation seeded NumPy globally, printed configured data layers,
and loaded PyVista for numerical work. Importing visualization changed PyVista's
global theme, and importing plot configuration replaced Matplotlib font settings.
Console initialization installed a process-wide display hook and recorded all
output indefinitely. Repeated logger lookups added duplicate handlers.

Zonation now imports rendering code only when requested. Random patterns use a
local or caller-supplied generator. Each plotter receives a local theme. Fresh
subprocess tests import from an arbitrary directory and verify quiet imports,
preserved random state, preserved global theme, and numerical imports without
PyVista.
Matplotlib styling now uses an explicit `plot_style()` context. Console output
is not recorded by default, the display hook is preserved, and logger lookups
reuse their existing handler. Import tests and output-capture tests cover these
behaviors.

### Medium: random pattern normalization sampled different values

The old pattern normalizer sampled a random field to calculate the scaling
factor, then sampled again to produce the output. The resulting weighted mean
did not match the reference. It also assumed a single cell block.

Patterns are now sampled once, normalized over all cells, and split back into
their blocks. A supplied `rng` makes random generation reproducible. Nonfinite
or zero weighted totals raise an error. Tests verify reproducibility and the
volume-weighted reference mean of 0.5. Exponential decrease now respects arbitrary
requested endpoints instead of assuming their sum is one.

### Medium: rendering and image APIs had surprising behavior

Iterator inputs were exhausted after the first frame, the time-course
`window_size` argument was ignored, and `conversion_factor` was unused. Plotters
were not reliably closed on setup errors. Image annotation required an external
`arial.ttf` and left its input image open. `vtu_to_xdmf` returned a mesh without
writing the requested conversion.

Rendering now retains iterator inputs, applies window settings and unit
conversion, and closes plotters in a `finally` block. Converted values use a copy
so the source mesh is preserved; derived color ranges use display units.
Annotation uses Pillow's portable default font or an explicit font path, creates
parents, and closes images. VTU conversion now writes readable XDMF/HDF output
with an explicit overwrite option. Tests cover iterator reuse, conversion values,
range ordering, plotter cleanup, annotation, and a VTU/XDMF roundtrip.

### Medium: HDF files remained open after failed XML serialization

Meshio's writer closes its HDF file after writing XML. If XML serialization
raises, that close is bypassed. A reproducer targeted a directory instead of an
XML file and confirmed that the HDF handle remained valid after the exception.

`XDMFTimeSeriesWriter` now closes the sidecar in `finally`, including this failure
path. The regression verifies the sidecar can immediately be reopened in truncate
mode after the metadata write fails.
GeoJSON plotting also closes its newly created figure if geometry parsing or
plot setup raises, without closing figures owned by the application.

### Medium: image montages held every input file and decoded image open

`merge_images` entered an image context for every panel and retained all of them
until saving the complete montage. A reproducer merging 100 one-pixel images
under a 32-descriptor process limit failed with `Too many open files`.

The helper now reads dimensions in a first pass, then opens, pastes, and closes
one image at a time. It retains the final canvas and current panel instead of
all decoded inputs and file handles. The same constrained workflow succeeds,
and a POSIX regression runs the file limit in a child process so other tests
and application state remain unaffected.

### Medium: directory output names could collide or vary by platform

Directory names were flattened by replacing `/` in a platform-specific path
string. Windows separators survived, and `a/b` could collide with `a__b`, silently
skipping or replacing one simulation. Every input path was also probed as if it
were a directory.

Discovery now visits VTK matches, deduplicates their parent directories, constructs
names from path components, and rejects collisions before conversion. Files in
the root input directory use that directory's name. Empty and missing sources
raise clear errors. A regression checks collision rejection without partial
output.

### Medium: research code and requirements shipped with the library

SPT, SimLivA, droplets, project geometry reconstruction, and workstation examples
lived under `porous_media.analyses`, so they shipped inside the library. Package
constants also inferred private data and result locations from the installation
directory. Spreadsheet and Arrow tools were direct runtime requirements even
though reusable library code did not use them.

All project workflows now live under repository-only `research/`, with their
synthetic tests under `research/tests/`. Historical mesh fixtures move unchanged.
Research is absent from wheels and source distributions, and library imports do
not depend on it. Research dependencies have a separate uv group, included by
`dev` for full-checkout checks. The library exposes bundled resources and accepts
explicit user paths. Packaging tests check the installed distribution boundary.

## Performance measurement

The same synthetic line mesh contains 2,000 cells: 500 periportal, 500 perivenous,
and 1,000 interior cells. Three calls use the default copy behavior on Python
3.14.4. The baseline median is **4.8172 seconds**; the spatial-index implementation
measures **0.00236 seconds**, approximately **2,000 times faster in this case**.
This compares the zonation operation, excluding imports and rendering. Geometry,
boundary distribution, hardware, and workload affect the result.

Reproduce the updated measurement:

```bash
uv run --no-sync python scripts/benchmark_zonation.py --cells 2000 --repeat 3
```

The implementation uses [SciPy KDTree](https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.KDTree.html)
for boundary queries. It avoids both the Python loop over every boundary for
each interior cell and an all-pairs distance matrix. Metadata and limit scans
also avoid `np.dstack` on each scalar element. No timing threshold is imposed on
CI because wall times vary across runners.

## Compatibility and operational limits

Migration instructions live in [research setup](../research.md),
[FEBio rendering](../febio.md), [zonation](../zonation.md), and
[data workflows](../xdmf.md). Notable behavior changes are explicit overwrite for
FEBio/VTU output, strict template parameters, XML escaping, safe export filenames,
real VTU conversion, removal of repository constants and the old analysis import
namespace, actual array shapes in metadata, and local random generators.
Matplotlib helper imports also no longer install styling globally; use the
explicit `plot_style()` context around figure creation.

Numerical NetCDF and JSON caches still require `overwrite=True` after external
source edits. They are local acceleration files, not independently versioned
datasets. Time series must retain fixed geometry and compatible fields. Native
VTK/HDF/XML readers remain external dependencies; this review does not make
arbitrary files or templates safe for an untrusted hosting service.

Historical research algorithms and private study selections remain project
code. This round separates and checks them without promoting their assumptions
into the public API or changing historical datasets and simulation results.

## Dependency advisory check

On 2026-10-19, the PyPA dependency auditor checked 42 applicable locked runtime
dependencies against PyPI's advisory service. It reported no known vulnerabilities
and skipped no applicable packages. The audit used a requirements export from
`uv.lock`, without resolving or installing a different dependency set. This is an
advisory-database result, separate from the source-level findings above.

```bash
uv export --locked --no-dev --no-emit-project --no-hashes \
  --format requirements-txt --output-file /tmp/porous-media-runtime-requirements.txt
uvx pip-audit --requirement /tmp/porous-media-runtime-requirements.txt \
  --no-deps --disable-pip --strict --progress-spinner off
```

## Validation

The final local checks pass: Ruff lint and formatting, ty with warnings treated
as errors, strict Zensical documentation, wheel and source builds, strict Twine
metadata checks, and all pre-commit hooks. Actionlint validates the workflow.

Installed-wheel tests pass all 99 library cases on Linux with Python 3.13.16 and
3.14.8. The four synthetic research tests pass independently. Package manifests
confirm that both release formats exclude research code while retaining bundled
HDF examples. Git records the historical mesh fixture moves with identical file
contents. No private study data or solver is required by these checks.
