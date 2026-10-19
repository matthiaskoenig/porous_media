# Installation

## Requirements

Use Python 3.13 or 3.14. CI tests Python 3.13 and 3.14 on Linux and Python 3.14
on Windows and macOS. Runtime dependencies, including meshio, HDF5, NumPy, SciPy,
PyVista, and Pillow, are installed with the package.

Python 3.15 is not yet supported. Its CI job is temporarily disabled while
upstream dependency compatibility remains incomplete. Use Python 3.13 or 3.14.
See [manual compatibility testing](development.md#required-checks) and
[issue #2](https://github.com/matthiaskoenig/porous_media/issues/2).

## Install a release

Create and activate a virtual environment, then install from PyPI:

```bash
python -m venv .venv
```

On Linux and macOS:

```bash
source .venv/bin/activate
```

On Windows PowerShell:

```powershell
.venv\Scripts\Activate.ps1
```

```bash
python -m pip install --upgrade pip
python -m pip install porous_media
python -c "import porous_media; print(porous_media.__version__)"
```

## Install the current source

To install the code on `main` without a development checkout:

```bash
python -m pip install git+https://github.com/matthiaskoenig/porous_media.git@main
```

For an editable checkout with the committed lockfile and development tools,
follow [development setup](development.md#environment).

## Optional external programs

### FFmpeg

Video and GIF helpers invoke `ffmpeg` from `PATH`. Install it through your
operating system package manager or the [FFmpeg download page](https://ffmpeg.org/download.html).
For example, on Debian or Ubuntu:

```bash
sudo apt install ffmpeg
ffmpeg -version
```

FFmpeg is not needed for XDMF processing, zonation export, or still images.

### Rendering

PyVista uses VTK and an OpenGL rendering backend. Interactive rendering needs
a display. Image export uses an off-screen plotter by default but still needs
a working VTK graphics backend. On servers, use a VTK-compatible EGL, OSMesa,
or virtual display setup for your environment.
See [PyVista installation](https://docs.pyvista.org/getting-started/installation.html)
for platform-specific details.

Start with the numerical [quickstart](quickstart.md) to verify installation
without opening a window. The [zonation example](zonation.md#bundled-example)
also supports `visualize=False`.

### FEBio

FEBio and FEBio Studio are installed independently. Python helpers render input
templates and process exported VTK files; they do not install FEBio, RoadRunner,
or research plugins. See [FEBio workflows](febio.md).
