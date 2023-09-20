.. image:: https://github.com/matthiaskoenig/porous_media/raw/main/docs/images/favicon/porous_media-100x100-300dpi.png
   :align: left
   :alt: porous_media logo

porous_media: python utilities for porous media analysis and visualization
==========================================================================

.. image:: https://github.com/matthiaskoenig/porous_media/workflows/CI-CD/badge.svg
   :target: https://github.com/matthiaskoenig/porous_media/workflows/CI-CD
   :alt: GitHub Actions CI/CD Status

.. image:: https://img.shields.io/pypi/v/porous_media.svg
   :target: https://pypi.org/project/porous_media/
   :alt: Current PyPI Version

.. image:: https://img.shields.io/pypi/pyversions/porous_media.svg
   :target: https://pypi.org/project/porous_media/
   :alt: Supported Python Versions

.. image:: https://img.shields.io/pypi/l/porous_media.svg
   :target: http://opensource.org/licenses/LGPL-3.0
   :alt: GNU Lesser General Public License 3

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.8335487.svg
   :target: https://doi.org/10.5281/zenodo.8335487
   :alt: Zenodo DOI

.. image:: https://img.shields.io/badge/code%20style-black-000000.svg
   :target: https://github.com/ambv/black
   :alt: Black

.. image:: http://www.mypy-lang.org/static/mypy_badge.svg
   :target: http://mypy-lang.org/
   :alt: mypy

porous_media is a collection of python utilities for working with porous media simulation results and meshes with source code available from 
`https://github.com/matthiaskoenig/porous_media <https://github.com/matthiaskoenig/porous_media>`__.

Features include among others

- Visualization of VTK results
- Image processing
- Mesh manipulation
- FEBio helpers
 
If you have any questions or issues please `open an issue <https://github.com/matthiaskoenig/porous_media/issues>`__.

How to cite
===========

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.8335487.svg
   :target: https://doi.org/10.5281/zenodo.8335487
   :alt: Zenodo DOI

Installation
============
`porous_media` is available from `pypi <https://pypi.python.org/pypi/porous_media>`__ and 
can be installed via:: 

    pip install porous_media

Required dependencies
----------------------
`porous_media` depends on `ffmpeg <https://ffmpeg.org/>`__ for creation of the videos 
which can be installed via on linux via::

    sudo apt-get install ffmpeg
    
Make sure you can execute the following successfully:: 

    ffmpeg

Develop version
---------------
The latest develop version can be installed via::

    pip install git+https://github.com/matthiaskoenig/porous_media.git@develop

Or via cloning the repository and installing via::

    git clone https://github.com/matthiaskoenig/porous_media.git
    cd porous_media
    pip install -e .

To install for development use::

    pip install -e .[development]


Contributing
============

Contributions are always welcome! Please read the `contributing guidelines
<https://github.com/matthiaskoenig/porous_media/blob/develop/.github/CONTRIBUTING.rst>`__ to
get started.

License
=======

* Source Code: `LGPLv3 <http://opensource.org/licenses/LGPL-3.0>`__
* Documentation: `CC BY-SA 4.0 <http://creativecommons.org/licenses/by-sa/4.0/>`__

The porous_media source is released under both the GPL and LGPL licenses version 2 or
later. You may choose which license you choose to use the software under.

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License or the GNU Lesser General Public
License as published by the Free Software Foundation, either version 2 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

Funding
=======
Matthias König is supported by the Federal Ministry of Education and Research (BMBF, Germany)
within the research network Systems Medicine of the Liver (**LiSyM**, grant number 031L0054) 
and by the German Research Foundation (DFG) within the Research Unit Programme FOR 5151 
"`QuaLiPerF <https://qualiperf.de>`__ (Quantifying Liver Perfusion-Function Relationship in Complex Resection - 
A Systems Medicine Approach)" by grant number 436883643 and by grant number 465194077 (Priority Programme SPP 2311, Subproject SimLivA).

© 2023 Matthias König
