![porous_media logo](https://raw.githubusercontent.com/matthiaskoenig/porous_media/main/docs/images/favicon/porous_media-100x100-300dpi.png)

# porous_media

[![CI](https://github.com/matthiaskoenig/porous_media/actions/workflows/main.yml/badge.svg)](https://github.com/matthiaskoenig/porous_media/actions/workflows/main.yml)
[![Version](https://img.shields.io/pypi/v/porous_media.svg)](https://pypi.org/project/porous_media/)
[![Python versions](https://img.shields.io/pypi/pyversions/porous_media.svg)](https://pypi.org/project/porous_media/)
[![MIT license](https://img.shields.io/pypi/l/porous_media.svg)](mit.txt)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8335487.svg)](https://doi.org/10.5281/zenodo.8335487)

Python utilities for porous media simulation results and meshes:

- Convert FEBio-exported VTK time series to XDMF/HDF5.
- Inspect and interpolate simulation data.
- Generate zonation patterns on annotated liver meshes.
- Render mesh data with PyVista and combine images and videos.
- Render FEBio input templates for parameter studies.

Python 3.13 or later is required. CI tests Python 3.13 and 3.14 on Linux and Python 3.14 on Windows and macOS.

## Installation

```bash
python -m pip install porous_media
```

Install the current source from `main`:

```bash
python -m pip install git+https://github.com/matthiaskoenig/porous_media.git@main
```

Video creation additionally needs [FFmpeg](https://ffmpeg.org/) on `PATH`.
Mesh conversion and numerical processing do not need FFmpeg or a display.
See the [installation guide](docs/installation.md) for rendering requirements.

## Quickstart

Save this script as `quickstart.py` and run `python quickstart.py`:

```python
from pathlib import Path

from porous_media import RESOURCES_DIR
from porous_media.data.xdmf_tools import XDMFInfo, vtks_to_xdmf

output = Path("results") / "example.xdmf"
vtks_to_xdmf(
    vtk_dir=RESOURCES_DIR / "vtk" / "vtk_timecourse",
    xdmf_path=output,
)
info = XDMFInfo.from_path(output)
print(f"{info.num_steps} steps: {info.tstart:g} to {info.tend:g}")
assert info.num_steps == 3
```

This uses bundled data and writes `example.xdmf`, `example.h5`, and
`example_limits.json` into `results/`. The summary prints `3 steps: 0 to 220`.

## Documentation and development

Start with the [documentation overview](docs/index.md), [user guides](docs/quickstart.md),
and [API reference](docs/api/index.md).
Documentation uses Zensical and API documentation is rendered from source docstrings.

```bash
git clone https://github.com/matthiaskoenig/porous_media.git
cd porous_media
uv sync --locked --group dev --group docs
uv run --no-sync pre-commit install
uv run --no-sync ruff check .
uv run --no-sync ruff format --check .
uv run --no-sync ty check
uv run --no-sync pytest
uv run --no-sync zensical build --clean --strict
```

See [CONTRIBUTING.md](CONTRIBUTING.md) for setup and pull requests,
[development](docs/development.md) for quality checks and repository policies,
and [RELEASE.md](RELEASE.md) for releases.
Project-specific scripts under `porous_media.analyses` may require unpublished data
and are not general-purpose examples.

## Citation and license

Cite the project using [doi:10.5281/zenodo.8335487](https://doi.org/10.5281/zenodo.8335487).
Source code is licensed under [MIT](mit.txt).
Documentation is licensed under [CC BY-SA 4.0](cc-by-sa-4.0.txt).

## Funding

Matthias König is supported by the German Research Foundation (DFG) within the Research Unit Programme FOR 5151
"QuaLiPerF (Quantifying Liver Perfusion-Function Relationship in Complex Resection -
A Systems Medicine Approach)" by grant number 436883643 and by grant number
465194077 (Priority Programme SPP 2311, Subproject SimLivA).

© 2023-2026 Matthias König
