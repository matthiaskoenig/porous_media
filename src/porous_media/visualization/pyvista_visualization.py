"""Visualization with pyvista.

Figure out how to add multiple layers to the visualization.
E.g. streamlines in addition to other data.

"""

from collections.abc import Generator, Iterable
from contextlib import contextmanager
from dataclasses import dataclass, field
from enum import Enum
from math import isfinite
from pathlib import Path

import meshio
import pyvista as pv
from pyvista.plotting._typing import ColormapOptions
from rich.progress import track

from porous_media._paths import filename_component
from porous_media.console import console
from porous_media.data.xdmf_tools import (
    AttributeType,
    DataLimits,
)
from porous_media.visualization.image_manipulation import merge_images


class VizType(str, Enum):  # noqa: UP042 - Preserve the public Enum string representation.
    """Type of visualization layer."""

    SCALAR = "Scalar"
    VECTOR_FIELD = "Vector"
    TENSOR = "Tensor"
    MATRIX = "Matrix"
    OTHER = "Other"


def _plot_theme() -> pv.themes.Theme:
    """Create a local theme without changing the application's PyVista settings."""
    theme = pv.themes.Theme()
    theme.window_size = [1000, 1000]
    theme.background = "white"
    theme.transparent_background = True
    theme.colorbar_orientation = "vertical"
    theme.font.family = "arial"
    theme.font.size = 20
    theme.font.title_size = 40
    theme.font.label_size = 30
    theme.font.color = "black"
    return theme


class DataRangeType(str, Enum):  # noqa: UP042 - Preserve the public Enum string representation.
    """Enum for handling the data range.

    Local: use the data range of every individual variable
    GLOBAL: use the global range over multiple variables
    """

    LOCAL = "LOCAL"
    GLOBAL = "GLOBAL"


@dataclass
class DataLayer:
    """Data layer for plot.

    Currently only a single data layer is supported but this will be extended
    to map multiple layers in a single visualiztion

    Attributes:
        sid: Name of the mesh data array to render.
        title: Scalar bar title.
        colormap: PyVista colormap name or definition.
        viz_type: Type of mesh attribute to activate.
        color_limits: Fixed minimum and maximum scalar values.
        scalar_bar: Whether to show a scalar bar.
        conversion_factor: Factor for converting data before plotting.
    """

    sid: str
    title: str
    colormap: ColormapOptions = "magma"
    viz_type: AttributeType = AttributeType.SCALAR
    color_limits: tuple[float, float] | None = None
    scalar_bar: bool = True
    conversion_factor: float | None = None

    def update_color_limits(
        self, data_limits: DataLimits, only_empty: bool = True
    ) -> None:
        """Update color limits.

        Args:
            data_limits: Available limits indexed by the data layer identifier.
            only_empty: Update only layers without existing limits when true.
        """
        if not only_empty or (only_empty and not self.color_limits):
            limits = data_limits.limits[self.sid]
            if self.conversion_factor is not None:
                if not isfinite(self.conversion_factor):
                    raise ValueError("conversion_factor must be finite")
                values = tuple(value * self.conversion_factor for value in limits)
                limits = (min(values), max(values))
            self.color_limits = limits


def _mesh_for_layer(mesh: pv.UnstructuredGrid, layer: DataLayer) -> pv.UnstructuredGrid:
    """Convert one layer without mutating the caller's mesh or other layers."""
    if layer.conversion_factor is None:
        return mesh
    if not isfinite(layer.conversion_factor):
        raise ValueError("conversion_factor must be finite")
    converted = mesh.copy(deep=True)
    converted[layer.sid] = converted[layer.sid] * layer.conversion_factor
    return converted


def visualize_datalayers_timecourse(
    xdmf_path: Path,
    output_dir: Path,
    data_layers: Iterable[DataLayer],
    window_size: tuple[int, int] = (600, 600),
) -> None:
    """Create visualizations for individual panels.

    Args:
        xdmf_path: Input XDMF time course.
        output_dir: Directory for rendered panels.
        data_layers: Layers to render separately for each time point.
        window_size: Width and height of each panel in pixels.
    """
    data_layers = tuple(data_layers)
    if not data_layers:
        raise ValueError("data_layers must contain at least one layer")
    settings = VisualizationSettings(window_size=list(window_size))

    # FIXME: Create more flexible visualization combining multiple datalayers into a single plot;
    # e.g., scalar value, vector field, streamlines, cloud visualization of steatosis.

    # create output dir
    if not output_dir.exists():
        output_dir.mkdir(parents=True)
        console.print(f"output_dir created: {output_dir}")

    with meshio.xdmf.TimeSeriesReader(xdmf_path) as reader:
        points, cells = reader.read_points_cells()

        tnum = reader.num_steps
        for k in track(
            range(tnum), description=f"Creating {tnum} panels for {xdmf_path.stem} ..."
        ):
            _, point_data, cell_data = reader.read_data(k)

            # Create mesh with single data point
            mesh: meshio.Mesh = meshio.Mesh(
                points=points, cells=cells, cell_data=cell_data, point_data=point_data
            )
            visualize_data_layers(
                mesh=mesh,
                data_layers=data_layers,
                image_name=f"sim_{k:05d}",
                output_dir=output_dir / "panels",
                visualization_settings=settings,
            )


@dataclass
class VisualizationSettings:
    """General visualization settings for panel.

    These should be determined once for a given geometry and then be used
    consistently.
    """

    # plotter settings
    window_size: list[int] = field(default_factory=lambda: [1000, 1000])
    off_screen: bool = (
        True  # False: blocking interactive visualization, True: batch mode
    )

    camera_position: tuple[float, float, float] = (0, 3e-4, 1e-3)
    zoom: float = 1.1


@contextmanager
def _plotter(
    settings: VisualizationSettings, *, interactive: bool = False
) -> Generator[pv.Plotter]:
    """Close each plotter even when mesh setup or screenshot creation fails."""
    plotter = pv.Plotter(
        window_size=settings.window_size,
        off_screen=False if interactive else settings.off_screen,
        theme=_plot_theme(),
    )
    try:
        yield plotter
    finally:
        plotter.close()


def visualize_data_layers(
    mesh: meshio.Mesh,
    data_layers: Iterable[DataLayer],
    output_dir: Path,
    image_name: str,
    visualization_settings: VisualizationSettings | None = None,
) -> None:
    """Visualize geometry with pyvista.

    Args:
        mesh: Geometry and data for one time point.
        data_layers: Layers to render as separate panels.
        output_dir: Directory for per-layer output folders.
        image_name: Output image name without its extension.
        visualization_settings: Camera and rendering settings, or defaults.
    """
    filename_component(image_name)
    data_layers = tuple(data_layers)
    if not data_layers:
        raise ValueError("data_layers must contain at least one layer")
    for layer in data_layers:
        filename_component(layer.sid)
        directory = output_dir / layer.sid
        if directory.is_symlink() or (directory / f"{image_name}.png").is_symlink():
            raise ValueError("Image output must not follow symbolic links")

    if not visualization_settings:
        # create default settings
        visualization_settings = VisualizationSettings()

    pvmesh = pv.utilities.from_meshio(mesh)

    # deactivate active sets
    pvmesh.set_active_tensors(None)
    pvmesh.set_active_scalars(None)
    pvmesh.set_active_vectors(None)

    # visualize data_layers
    data_layers_dict: dict[str, DataLayer] = {dl.sid: dl for dl in data_layers}
    for name, data_layer in data_layers_dict.items():
        layer_mesh = _mesh_for_layer(pvmesh, data_layer)
        with _plotter(visualization_settings) as p:
            if data_layer.viz_type == "Scalar":
                layer_mesh.set_active_scalars(name=name)
            elif data_layer.viz_type == "Vector":
                layer_mesh.set_active_vectors(name=name)
            elif data_layer.viz_type == "Tensor":
                layer_mesh.set_active_tensors(name=name)

            # FIXME: add multiple layers on top of each other
            actor = p.add_mesh(
                layer_mesh,
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
                # PyVista's functools.wraps exposes an unbound method signature to ty.
                p.add_scalar_bar(  # ty: ignore[missing-argument]
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
    ncols: int | None = None,
    nrows: int | None = None,
) -> list[Path]:
    """Create combined images for all timepoints."""
    filename_component(direction)
    selection = tuple(filename_component(name) for name in selection)
    if not selection:
        raise ValueError("selection must contain at least one layer")
    if direction not in {"vertical", "horizontal", "square", "custom"}:
        raise ValueError(f"Unsupported image direction: {direction!r}")
    if num_steps < 0:
        raise ValueError("num_steps must be nonnegative")
    image_dir: Path = output_dir / direction
    image_dir.mkdir(parents=True, exist_ok=True)
    all_images: list[Path] = []
    for k in track(range(num_steps), description="Creating combined images ..."):
        images: list[Path] = []
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
    visualization_settings: VisualizationSettings,
) -> None:
    """Visualize geometry with pyvista."""
    pvmesh = pv.utilities.from_meshio(mesh)
    pvmesh = _mesh_for_layer(pvmesh, data_layer)

    # name of variable
    name = data_layer.sid

    # deactivate active sets
    pvmesh.set_active_tensors(None)
    pvmesh.set_active_scalars(None)
    pvmesh.set_active_vectors(None)

    # visualize data_layer
    with _plotter(visualization_settings, interactive=True) as p:
        if data_layer.viz_type == "Scalar":
            pvmesh.set_active_scalars(name=name)
        elif data_layer.viz_type == "Vector":
            pvmesh.set_active_vectors(name=name)
        elif data_layer.viz_type == "Tensor":
            pvmesh.set_active_tensors(name=name)

        # PyVista's functools.wraps exposes an unbound method signature to ty.
        p.show_grid()  # ty: ignore[missing-argument]

        # FIXME: add multiple layers on top of each other, combination of data layers
        actor = p.add_mesh(
            pvmesh,
            show_vertices=True,
            show_edges=True,
            render_points_as_spheres=True,
            point_size=3,
            line_width=1.0,
            cmap=data_layer.colormap,
            show_scalar_bar=False,
            edge_color="darkgray",
            # specular=0.5, specular_power=15,
            clim=data_layer.color_limits,
            # pbr=True,
            # metallic=0.8,
            # roughness=0.0,
            # diffuse=1
        )
        # set up lighting
        # light = pv.Light((-2, 2, 0), (0, 0, 0), 'white')
        # p.add_light(light)
        #
        light = pv.Light((-2, 0, 0), (0, 0, 0), (0.7, 0.0862, 0.0549))
        p.add_light(light)
        #
        # light = pv.Light((0, 0, 10), (0, 0, 0), 'white')
        # p.add_light(light)

        # p.add_floor('-z', lighting=True, color='white', pad=1.0)
        p.add_floor(  # ty: ignore[missing-argument]
            face="-z", lighting=True, color="white", pad=0.5
        )
        p.enable_shadows()  # ty: ignore[missing-argument]

        # TODO:
        # # arrows for vectors: https://docs.pyvista.org/version/stable/examples/01-filter/glyphs.html
        # p.add_mesh(
        #     # vector on point data
        #     pvmesh.arrows,
        #     lighting=False,
        # )

        if data_layer.scalar_bar:
            p.add_scalar_bar(  # ty: ignore[missing-argument]
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


def xdmf_to_mesh(xdmf_path: Path, k: int = 0) -> meshio.Mesh:
    """XDMF to mesh."""
    with meshio.xdmf.TimeSeriesReader(xdmf_path) as reader:
        points, cells = reader.read_points_cells()
        _, point_data, cell_data = reader.read_data(k)

        return meshio.Mesh(
            points=points, cells=cells, cell_data=cell_data, point_data=point_data
        )
