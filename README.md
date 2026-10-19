![porous_media logo](https://raw.githubusercontent.com/matthiaskoenig/porous_media/main/docs/images/favicon/porous_media-100x100-300dpi.png)

# porous_media: Python tools for porous media

[![CI-CD](https://github.com/matthiaskoenig/porous_media/actions/workflows/ci-cd.yml/badge.svg)](https://github.com/matthiaskoenig/porous_media/actions/workflows/ci-cd.yml)
[![Documentation](https://img.shields.io/badge/docs-porous_media-3f51b5.svg)](https://matthiaskoenig.github.io/porous_media/)
[![Version](https://img.shields.io/pypi/v/porous_media.svg)](https://pypi.org/project/porous_media/)
[![Python versions](https://img.shields.io/pypi/pyversions/porous_media.svg)](https://pypi.org/project/porous_media/)
[![MIT License](https://img.shields.io/pypi/l/porous_media.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8335487.svg)](https://doi.org/10.5281/zenodo.8335487)

`porous_media` provides Python utilities for porous media simulation results and meshes:

- Convert FEBio-exported VTK time series to XDMF/HDF5.
- Inspect and interpolate simulation data.
- Generate zonation patterns on annotated liver meshes.
- Render mesh data with PyVista and combine images and videos.
- Render FEBio input templates for parameter studies.

The documentation is available at **[matthiaskoenig.github.io/porous_media](https://matthiaskoenig.github.io/porous_media/)**:

- [Installation](https://matthiaskoenig.github.io/porous_media/installation/) and [quickstart](https://matthiaskoenig.github.io/porous_media/quickstart/)
- [User guides](https://matthiaskoenig.github.io/porous_media/xdmf/) and [API reference](https://matthiaskoenig.github.io/porous_media/api/)
- [Contributing](https://matthiaskoenig.github.io/porous_media/contributing/) and [development](https://matthiaskoenig.github.io/porous_media/development/)
- [Release guide](https://matthiaskoenig.github.io/porous_media/release/) and [release notes](https://matthiaskoenig.github.io/porous_media/release-notes/)

Report bugs or request features through the [issue tracker](https://github.com/matthiaskoenig/porous_media/issues).

## How to cite

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8335487.svg)](https://doi.org/10.5281/zenodo.8335487)

If you use `porous_media` please cite the archived software on [Zenodo](https://doi.org/10.5281/zenodo.8335487):

> König, M. (2026). *porous_media: python utilities for porous media analysis and visualization* (Version 0.4.0) [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.23471452

```bibtex
@software{konig_porous_media,
  author    = {König, Matthias},
  title     = {porous_media: python utilities for porous media analysis and visualization},
  year      = {2026},
  month     = oct,
  version   = {0.4.0},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.23471452},
  url       = {https://doi.org/10.5281/zenodo.23471452},
}
```

Citation metadata is available in [CITATION.cff](CITATION.cff).

## License

- Source Code: [MIT](https://opensource.org/license/MIT)
- Documentation: [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/)

## Funding

Matthias König is supported by the German Research Foundation (DFG) within the Research Unit Programme FOR 5151
"QuaLiPerF (Quantifying Liver Perfusion-Function Relationship in Complex Resection -
A Systems Medicine Approach)" by grant number 436883643 and by grant number
465194077 (Priority Programme SPP 2311, Subproject SimLivA).

© 2023-2026 Matthias König, [https://livermetabolism.com](https://livermetabolism.com)
