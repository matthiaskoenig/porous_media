# Visualization and video

PyVista renders mesh fields into individual image panels. Each `DataLayer`
selects a field by its exact mesh identifier. A layer produces a separate panel;
the helper does not overlay multiple fields in one image.

## Render a static mesh

This example uses the bundled zonation mesh and needs a working VTK rendering
backend. Run the numerical [quickstart](quickstart.md) first when checking a
new installation.

```python
from pathlib import Path

from porous_media import RESOURCES_DIR
from porous_media.mesh.mesh_zonation import ZonatedMesh, ZonationPatterns
from porous_media.visualization.pyvista_visualization import (
    DataLayer,
    VisualizationSettings,
    visualize_data_layers,
)

mesh = ZonatedMesh.create_zonated_mesh_from_vtk(
    vtk_path=RESOURCES_DIR / "zonation" / "mesh_zonation_lobulus.vtk",
    patterns=[ZonationPatterns.linear_increase],
)
visualize_data_layers(
    mesh=mesh,
    data_layers=[DataLayer(sid="position", title="Periportal to perivenous")],
    output_dir=Path("results") / "images",
    image_name="zonation",
    visualization_settings=VisualizationSettings(off_screen=True),
)
```

The image is written to `results/images/position/zonation.png`.
The default camera fits the bundled liver geometry. Adjust `camera_position`
and `zoom` in `VisualizationSettings` for other mesh scales and orientations.
Off-screen mode still requires a graphics backend; it does not remove VTK's
OpenGL requirement.

## Render a time series

After the [quickstart](quickstart.md), render its scalar pressure field:

```python
from pathlib import Path

from porous_media.data.xdmf_tools import DataLimits
from porous_media.visualization.pyvista_visualization import (
    DataLayer,
    visualize_datalayers_timecourse,
)

xdmf_path = Path("results") / "quickstart" / "example.xdmf"
layer = DataLayer(
    sid="effective_fluid_pressure_TPM",
    title="Effective fluid pressure",
    color_limits=DataLimits.from_xdmf(xdmf_path).limits["effective_fluid_pressure_TPM"],
)
visualize_datalayers_timecourse(
    xdmf_path=xdmf_path,
    output_dir=Path("results") / "pressure",
    data_layers=[layer],
)
```

This creates `sim_00000.png`, `sim_00001.png`, and `sim_00002.png` under
`results/pressure/panels/effective_fluid_pressure_TPM/`.
Use fixed color limits to compare frames consistently. Confirm field units
from the simulation model; the rendering helper does not infer physical units.

## Combine images and create videos

`merge_images` combines existing panels vertically, horizontally, or in a grid.
`annotate_image_text` adds text to an existing image. See the
[visualization API](api/visualization.md) for their parameters.

With FFmpeg installed, combine the pressure frames into a video:

```python
from pathlib import Path

from porous_media.visualization.video import create_gif_from_video, create_video

output_dir = Path("results") / "pressure"
video = output_dir / "pressure.mp4"
create_video(
    image_pattern=str(
        output_dir / "panels" / "effective_fluid_pressure_TPM" / "sim_%05d.png"
    ),
    video_path=video,
    frame_rate=10,
)
create_gif_from_video(video_path=video, gif_path=video.with_suffix(".gif"))
```

The frame rate controls playback speed. It does not preserve irregular physical
time intervals between simulation steps. The video helpers overwrite outputs
by default. Confirm output files and FFmpeg diagnostics after encoding.

## Common problems

Batch exports use data-layer identifiers and image names as filename components.
Choose portable names without path separators, traversal, or Windows device names.
Unsafe paths are rejected before rendering. Plotters close even when setup or
screenshot creation fails, and importing the module preserves the application's
global PyVista theme. Iterator inputs are retained across time steps.

`DataLayer.conversion_factor` converts displayed values without changing the
source mesh. Ranges set through `update_color_limits` use the same display units.
Explicit color limits should already use those units.

Annotation uses Pillow's bundled default font. Pass `font_path` and `font_size`
to select an installed font explicitly. The helper closes its input image and
creates the output directory when needed.

Montage generation opens one input image at a time, so large panel collections
do not retain one file descriptor and decoded buffer per input. Optional
Matplotlib styling is explicit and restores prior settings:

```python
import matplotlib.pyplot as plt

from porous_media.visualization.plots import plot_style

with plot_style():
    figure, axes = plt.subplots()
    try:
        axes.plot([0, 1], [0, 1])
        figure.savefig("line.png")
    finally:
        plt.close(figure)
```

`porous_media.visualization.vtu_reader.vtu_to_xdmf` now writes the requested static
XDMF/HDF conversion and returns its mesh. This replaces its previous read-only
behavior. Outputs are preserved unless `overwrite=True` is supplied.

- **Field not found:** inspect `XDMFInfo.point_data`, `XDMFInfo.cell_data`, or
  `mesh.point_data` and `mesh.cell_data` before choosing a layer.
- **Blank or clipped mesh:** check camera settings and coordinate scale.
- **OpenGL or display errors:** configure a VTK-compatible graphics backend;
  see [installation](installation.md#rendering).
- **FFmpeg not found:** install FFmpeg and verify `ffmpeg -version` in the same
  environment used to launch Python.
