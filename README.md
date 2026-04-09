![porous_media logo](https://github.com/matthiaskoenig/porous_media/raw/main/docs/images/favicon/porous_media-100x100-300dpi.png)

# porous_media: python utilities for SBML
[![GitHub Actions CI/CD Status](https://github.com/matthiaskoenig/porous_media/workflows/CI-CD/badge.svg)](https://github.com/matthiaskoenig/porous_media/actions/workflows/main.yml)
[![Version](https://img.shields.io/pypi/v/porous_media.svg)](https://pypi.org/project/porous_media/)
[![Python Versions](https://img.shields.io/pypi/pyversions/porous_media.svg)](https://pypi.org/project/porous_media/)
[![MIT License](https://img.shields.io/pypi/l/porous_media.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8335487.svg)](https://doi.org/10.5281/zenodo.8335487)

porous_media is a collection of python utilities for working with porous media simulation results and meshes with source code available from
[https://github.com/matthiaskoenig/porous_media](https://github.com/matthiaskoenig/porous_media).

Features include among others

- Visualization of VTK results
- Image processing
- Mesh manipulation
- FEBio helpers

If you have any questions or issues please [open an issue](https://github.com/matthiaskoenig/porous_media/issues).


## How to cite
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8335487.svg)](https://doi.org/10.5281/zenodo.8335487)

## License
- Source Code: [MIT](https://opensource.org/license/MIT)
- Documentation: [CC BY-SA 4.0](http://creativecommons.org/licenses/by-sa/4.0/)

## Installation
### Dependencies
porous_media depends on [`ffmpeg`](https://ffmpeg.org/) for creation of the videos
which can be installed via on linux via
```bash
sudo apt -y install ffmpeg
```

Make sure you can execute the following successfully from the command line:
```bash
ffmpeg
```

### pip
porous_media is available from [pypi](https://pypi.python.org/pypi/porous_media) and
can be installed via
```bash
pip install porous_media
```

The develop version can be installed via
```bash
pip install git+https://github.com/matthiaskoenig/porous_media.git@develop
```

## Development
Setup develop environment
```bash
uv sync --group dev
```
Setup pre-commit
```bash
pre-commit install
pre-commit run
```

## Testing
Run all tests in parallel
```bash
tox run-parallel
```

Run single target
```bash
tox r -e py314
```


## Funding
Matthias König is supported and by the German Research Foundation (DFG) within the Research Unit Programme FOR 5151
"QuaLiPerF (Quantifying Liver Perfusion-Function Relationship in Complex Resection -
A Systems Medicine Approach)" by grant number 436883643 and by grant number
465194077 (Priority Programme SPP 2311, Subproject SimLivA).

© 2023-2026 Matthias König
