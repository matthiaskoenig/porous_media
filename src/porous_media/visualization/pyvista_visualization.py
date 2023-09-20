"""Visualization with pyvista.

Figure out how to add multiple layers to the visualization.
E.g. streamlines in addition to other data.

"""
import tempfile
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import meshio
import pyvista as pv
from rich.progress import track

from porous_media import RESOURCES_DIR, RESULTS_DIR
from porous_media.console import console
from porous_media.data.xdmf_tools import DataLimits, XDMFInfo, vtks_to_xdmf
from porous_media.mesh.mesh_tools import mesh_to_xdmf
from porous_media.visualization.image_manipulation import merge_images


# global configuration
pv.global_theme.window_size = [1000, 1000]
pv.global_theme.background = "white"
pv.global_theme.transparent_background = True
# pv.global_theme.cmap = 'RdBu'
pv.global_theme.colorbar_orientation = "vertical"

pv.global_theme.font.family = "arial"
pv.global_theme.font.size = 20
pv.global_theme.font.title_size = 40
pv.global_theme.font.label_size = 30
pv.global_theme.font.color = "black"


class DataRangeType(str, Enum):
    """Enum for handling the data range.

    Local: use the data range of every individual variable
    GLOBAL: use the global range over multiple variables
    """

    LOCAL = "LOCAL"
    GLOBAL = "GLOBAL"


@dataclass
class DataLayer:
    """Data layer for plot with all information to visualize the layer.

    This will be most likely extended to account for
    order of plotting, opacity and so on.

    :param conversion_factor: factor f to convert data before plotting with data_new = data * f
    """

    sid: str
    title: str
    data_type: str  # Scalar, Vector, Tensor
    colormap: str
    color_limits: Optional[Tuple[float, float]] = None
    scalar_bar: bool = True
    # plot_type: str
    conversion_factor: Optional[float] = None

    def update_color_limits(self, data_limits: DataLimits):
        self.color_limits = data_limits.limits[self.sid]


def visualize_datalayers_timecourse(
    xdmf_path: Path,
    output_dir: Path,
    data_layers: Iterable[DataLayer],
    window_size: Tuple[int, int] = (600, 600),
) -> None:
    """Create visualizations for individual panels.

    :param xdmf_path: timecourse xdmf
    :param data_layers: iterable of data layers to visualize, for every layer a plot is generated
    """

    # FIXME: Create more flexible visualization combining multiple datalayers into a single plot;
    # e.g., scalar value, vector field, streamlines, cloud visualization of steatosis.

    # create output dir
    if not output_dir.exists():
        output_dir.mkdir(parents=True)
        console.print(f"output_dir created: {output_dir}")

    with meshio.xdmf.TimeSeriesReader(xdmf_path) as reader:
        points, cells = reader.read_points_cells()

        tnum = reader.num_steps
        for k in track(range(tnum), description=f"Creating {tnum} panels for {xdmf_path.stem} ..."):
            t, point_data, cell_data = reader.read_data(k)

            # Create mesh with single data point
            mesh: meshio = meshio.Mesh(
                points=points, cells=cells, cell_data=cell_data, point_data=point_data
            )
            visualize_data_layers(
                mesh=mesh,
                data_layers=data_layers,
                image_name=f"sim_{k:05d}",
                output_dir=output_dir / "panels",
                window_size=window_size,
            )


@dataclass
class VisualizationSettings:
    """General visualization settings for panel.

    These should be determined once for a given geometry and then be used
    consistently.
    """
    # plotter settings
    window_size: Tuple[float, float] = (1000, 1000)
    off_screen: bool = True  # False: blocking interactive visualization, True: batch mode

    camera_position: Tuple[float, float, float] = (0, 3e-4, 1e-3)
    zoom: float = 1.1



def visualize_data_layers(
    mesh: meshio.Mesh,
    data_layers: Iterable[DataLayer],
    output_dir: Path,
    image_name: str,
    visualization_settings: Optional[VisualizationSettings]
) -> None:
    """Visualize geometry with pyvista.

    :param mesh: mesh with single time point scalar data
    :param image_name: name of the created image, without extension.
    """
    if not visualization_settings:
        # create default settings
        visualization_settings = VisualizationSettings()

    # FIXME: make this the function for a single plot

    # create VTK from mesh to read for pyvista
    # FIXME: better handling of grid
    xdmf_tmp = tempfile.NamedTemporaryFile(suffix=".xdmf")
    mesh_to_xdmf(m=mesh, xdmf_path=Path(xdmf_tmp.name), test_read=False)
    grid = pv.read(xdmf_tmp.name)

    # deactivate active sets
    grid.set_active_tensors(None)
    grid.set_active_scalars(None)
    grid.set_active_vectors(None)

    # visualize data_layers
    data_layers_dict: Dict[str, DataLayer] = {dl.sid: dl for dl in data_layers}
    for name, data_layer in data_layers_dict.items():
        p = pv.Plotter(
            window_size=visualization_settings.window_size,
            # title="TPM",
            off_screen=visualization_settings.off_screen,
        )
        if data_layer.data_type == "Scalar":
            grid.set_active_scalars(name=name)
        elif data_layer.data_type == "Vector":
            grid.set_active_vectors(name=name)
        elif data_layer.data_type == "Tensor":
            grid.set_active_tensors(name=name)

        # FIXME: add multiple layers on top of each other
        actor = p.add_mesh(
            grid,
            show_edges=True,
            render_points_as_spheres=True,
            point_size=3,
            show_vertices=True,
            line_width=1.0,
            cmap=data_layer.colormap,
            show_scalar_bar=False,
            edge_color="darkgray",
            # specular=0.5, specular_power=15,
            clim=data_layer.color_limits,
        )
        if data_layer.scalar_bar:
            p.add_scalar_bar(
                title=data_layer.title,
                n_labels=5,
                bold=True,
                # height=0.6,
                width=0.6,
                vertical=False,
                position_y=0.0,
                position_x=0.2,
                mapper=actor.mapper,
                fmt="%.1f",
            )

            # set the color limits
            if data_layer.color_limits is not None:
                p.update_scalar_bar_range(
                    clim=data_layer.color_limits,
                    name=data_layer.title,
                )

        # Camera position to zoom to face
        p.camera_position = visualization_settings.camera_position
        p.camera.zoom(visualization_settings.zoom)

        output_subdir = Path(output_dir) / f"{name}"
        output_subdir.mkdir(exist_ok=True, parents=True)

        # FIXME: have a look at the other serialization formats such as interactive HTML
        p.show(
            # cpos=top_view,
            screenshot=Path(output_subdir / f"{image_name}.png")
        )


def create_combined_images(
    num_steps: int,
    output_dir: Path,
    selection: Iterable[str],
    direction: str,
    ncols: Optional[int] = None,
    nrows: Optional[int] = None,
) -> List[Path]:
    """Create combined images for all timepoints."""

    image_dir: Path = output_dir / direction
    image_dir.mkdir(parents=True, exist_ok=True)
    all_images: List[Path] = []
    for k in track(range(num_steps), description="Creating combined images ..."):
        images: List[Path] = []
        for name in selection:
            img_path = output_dir / "panels" / name / f"sim_{k:05d}.png"
            images.append(img_path)

        combined_image: Path = output_dir / direction / f"sim_{k:05d}.png"
        merge_images(
            paths=images,
            direction=direction,
            output_path=combined_image,
            ncols=ncols,
            nrows=nrows,
        )

        all_images.append(combined_image)

    return all_images


def visualize_interactive(
    mesh: meshio.Mesh,
    data_layer: DataLayer,
    visualization_settings: VisualizationSettings
) -> None:
    """Visualize geometry with pyvista."""
    xdmf_tmp = tempfile.NamedTemporaryFile(suffix=".xdmf")
    mesh_to_xdmf(m=mesh, xdmf_path=Path(xdmf_tmp.name), test_read=False)
    grid = pv.read(xdmf_tmp.name)

    # name of variable
    name = data_layer.sid

    # deactivate active sets
    grid.set_active_tensors(None)
    grid.set_active_scalars(None)
    grid.set_active_vectors(None)

    # visualize data_layer
    p = pv.Plotter(
        window_size=visualization_settings.window_size,
    )
    if data_layer.data_type == "Scalar":
        grid.set_active_scalars(name=name)
    elif data_layer.data_type == "Vector":
        grid.set_active_vectors(name=name)
    elif data_layer.data_type == "Tensor":
        grid.set_active_tensors(name=name)

    # FIXME: add multiple layers on top of each other, combination of data layers
    actor = p.add_mesh(
        grid,
        show_edges=True,
        render_points_as_spheres=True,
        point_size=3,
        show_vertices=True,
        line_width=1.0,
        cmap=data_layer.colormap,
        show_scalar_bar=False,
        edge_color="darkgray",
        # specular=0.5, specular_power=15,
        clim=data_layer.color_limits,
    )

    # arrows for vectors: https://docs.pyvista.org/version/stable/examples/01-filter/glyphs.html
    p.add_mesh(
        # vector on point data
        grid.arrows,
        lighting=False,
    )

    if data_layer.scalar_bar:
        p.add_scalar_bar(
            title=data_layer.title,
            n_labels=5,
            bold=True,
            # height=0.6,
            width=0.6,
            vertical=False,
            position_y=0.0,
            position_x=0.2,
            mapper=actor.mapper,
            fmt="%.1f",
        )

        # set the color limits
        if data_layer.color_limits is not None:
            p.update_scalar_bar_range(
                clim=data_layer.color_limits,
                name=data_layer.title,
            )

    # Camera position to zoom to face
    p.camera_position = visualization_settings.camera_position
    p.camera.zoom(visualization_settings.zoom)

    p.show()


def xdmf_to_mesh(xdmf_path, k: int = 0) -> meshio.Mesh:
    """XDMF to mesh."""
    with meshio.xdmf.TimeSeriesReader(xdmf_path) as reader:
        points, cells = reader.read_points_cells()
        t, point_data, cell_data = reader.read_data(k)

        return meshio.Mesh(
            points=points, cells=cells, cell_data=cell_data, point_data=point_data
        )


if __name__ == "__main__":
    # FIXME: make this in an example and add path to resources

    # Simple test that all variables are read

    # vtk_dir = RESOURCES_DIR / "vtk" / "vtk_single"
    vtk_dir = RESOURCES_DIR / "vtk" / "vtk_timecourse"
    xdmf_path = RESULTS_DIR / "vtk_test.xdmf"
    vtks_to_xdmf(vtk_dir, xdmf_path=xdmf_path, overwrite=True)
    xdmf_info: XDMFInfo = XDMFInfo.from_path(xdmf_path)
    console.print(xdmf_info)

    data_layers: List[DataLayer] = [
        DataLayer(
            sid="displacement",
            title="displacement",
            colormap="magma",
            data_type="Vector",
        ),
        DataLayer(
            sid="fluid_flux_TPM",
            title="fluid flux [m/s]",
            colormap="magma",
            data_type="Vector",
        ),
        DataLayer(
            sid="pressure", title="Pressure", colormap="magma", data_type="Scalar"
        ),
        DataLayer(
            sid="effective_fluid_pressure_TPM",
            title="Effective fluid pressure [Pa]",
            colormap="magma",
            data_type="Scalar",
        ),
        DataLayer(
            sid="rr_(S)", title="Substrate S [mM]", colormap="magma", data_type="Scalar"
        ),
    ]
    # visualize_datalayers_timecourse(
    #     xdmf_path=xdmf_path,
    #     data_layers=data_layers,
    #     output_dir=RESULTS_DIR / "vtk_test_images"
    # )

    mesh = xdmf_to_mesh(xdmf_path, k=1)
    visualize_interactive(
        mesh,
        data_layer=data_layers[0],
    )
    # FIXME: support calculation of new variables for visualization
