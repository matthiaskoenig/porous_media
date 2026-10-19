# Research archive

Project workflows run from a repository checkout with their required datasets
and configuration. Use the [installation guide](installation.md),
[FEBio guide](febio.md), and [development guide](development.md) for supported
package workflows.

## Separate research subspace

The repository keeps reusable code and project workflows in different locations:

| Location | Responsibility | Published with the library |
| --- | --- | --- |
| `src/porous_media/` | Reusable data, mesh, visualization, and FEBio APIs | Yes |
| `src/porous_media/resources/` | Small bundled examples | Yes |
| `research/spt/`, `research/simliva/`, `research/droplets/` | Project workflows | No |
| `research/geometry_conversion.py` | Study-specific lobule reconstruction | No |
| `research/examples/` | Archived workstation scripts and examples | No |
| `research/tests/` | Synthetic research regressions | No |
| `data/`, `results/` | Local datasets and generated results | No |

Research imports the library, while the library never imports research. Both
subspaces share the committed uv lockfile and remain covered by Ruff and ty.
Research dependencies are declared in the `research` group, also included by
`dev` for full-checkout validation. A study-only checkout can use
`uv sync --locked --group research --no-dev`.

## Configure and run a study

Run project modules from the repository root after supplying their required
datasets and simulation metadata. Configuration is explicit:

- `POROUS_MEDIA_DATA_DIR` defaults to the checkout's `data/` directory.
- `POROUS_MEDIA_RESULTS_DIR` defaults to the checkout's `results/` directory.
- `POROUS_MEDIA_SPT_DATA` defaults to `DATA_DIR / "spt"`.

For example, after the development environment is synced:

```bash
export POROUS_MEDIA_SPT_DATA=/path/to/private/spt/data
uv run --no-sync python -m research.spt.information
```

On PowerShell, set `$env:POROUS_MEDIA_SPT_DATA` instead of using `export`.
Other historical scripts can retain study-specific selections and paths; review
those settings before running them. Importing `research` does not create folders
or run analyses.

## Import migration

Replace `porous_media.analyses.spt` imports with `research.spt` and
`porous_media.analyses.geometry_conversion` with `research.geometry_conversion`.
Use module execution from the checkout rather than running a nested file directly.
The old namespace is removed from the installed package; there is no compatibility
alias that pulls research back into library distributions.

The library exports `RESOURCES_DIR` for bundled inputs. Its former repository
constants `BASE_DIR`, `DATA_DIR`, `RESULTS_DIR`, and `EXAMPLE_*` are removed.
Pass explicit paths in library applications. Research scripts can import
`BASE_DIR`, `DATA_DIR`, and `RESULTS_DIR` from `research`.

Run `uv run --no-sync pytest -q research/tests` for portable synthetic coverage.
Neither these tests nor the library tests need private study datasets or FEBio.

## Study documentation

The [SPT guide](research/spt.md) describes the complete resumable workflow for
cataloguing, reconstructing, analyzing, and visualizing simulation and convergence
datasets. Run `uv run --no-sync python -m research.spt.pipeline` after supplying
the existing solver results. Numerical analysis defaults to homogeneous lobules;
mixed sectors are reserved for visualization. SPT entry points share explicit
configuration and serialized results rather than historical dataset selections.
