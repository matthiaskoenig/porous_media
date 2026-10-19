# porous_media

`porous_media` provides Python utilities for porous media simulation results and
meshes. It connects FEBio-exported VTK data with XDMF time series, mesh zonation,
and PyVista visualizations.

## Start here

1. [Install the package](installation.md) on Python 3.13 or later.
2. [Run the quickstart](quickstart.md) with bundled simulation data.
3. [Convert and inspect XDMF](xdmf.md), [generate zonation](zonation.md), or
   [render mesh data](visualization.md).
4. Use the [FEBio guide](febio.md) to prepare inputs and process solver outputs.

## Workflow

| Task | Input | Output |
| --- | --- | --- |
| Convert a time series | FEBio-exported legacy VTK files | XDMF, HDF5, and cached limits |
| Interpolate a time series | XDMF with a fixed mesh | XDMF at requested times |
| Generate zonation | Mesh with boundary labels and cell volumes | Cell positions and zonation patterns |
| Render results | Mesh data or an XDMF time series | Image panels |
| Create a video | Numbered images | Video or GIF using FFmpeg |
| Prepare a parameter study | Jinja2 FEBio template and contexts | Rendered `.feb` files |

The package does not run a finite element solver. FEBio, its plugins, and FFmpeg
are external programs. Numerical data processing does not require a display.

## Scope and reference

The [API reference](api/index.md) documents the reusable data, mesh,
visualization, and FEBio modules. Scripts in `porous_media.analyses` are research
workflows with project-specific paths and data requirements.
They are not a supported command-line interface or portable tutorials.

The [research archive](research.md) preserves historical setup and study notes separately from the current
user guides.

Development uses [uv](https://docs.astral.sh/uv/), Ruff, ty, pytest, and Zensical.
See [contributing](contributing.md) for bug reports and pull requests, and
[development](development.md) for the locked environment and required checks.

## Project

- [Source and issues](https://github.com/matthiaskoenig/porous_media)
- [Package on PyPI](https://pypi.org/project/porous_media/)
- [Citation: doi:10.5281/zenodo.8335487](https://doi.org/10.5281/zenodo.8335487)
- [MIT source license](https://github.com/matthiaskoenig/porous_media/blob/main/LICENSE)
- [CC BY-SA 4.0 documentation license](cc-by-sa-4.0.txt)
