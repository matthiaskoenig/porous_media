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

from porous_media.console import console
from porous_media.data.xdmf_tools import DataLimits
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
    """

    sid: str
    title: str
    data_type: str  # Scalar, Vector, Tensor
    colormap: str
    color_limits: Optional[Tuple[float, float]] = None
    scalar_bar: bool = True
    # plot_type: str

    def update_color_limits(self, data_limits: DataLimits):
        self.color_limits = data_limits.limits[self.sid]


def visualize_datalayers_timecourse(
    xdmf_path: Path,
    output_dir: Path,
    data_layers: List[DataLayer],
    window_size: Tuple[int, int] = (600, 600),
) -> None:
    """Create visualizations for individual panels."""
    # create output dir
    if not output_dir.exists():
        output_dir.mkdir(parents=True)
        console.print(f"output_dir created: {output_dir}")

    with meshio.xdmf.TimeSeriesReader(xdmf_path) as reader:
        points, cells = reader.read_points_cells()

        tnum = reader.num_steps
        for k in track(range(tnum), description="Create panels for data layers ..."):
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


def visualize_data_layers(
    mesh: meshio.Mesh,
    data_layers: List[DataLayer],
    output_dir: Path,
    image_name: str,
    window_size: Tuple[float, float] = (1000, 1000),
) -> None:
    """Visualize geometry with pyvista.

    :param mesh: mesh with single time point scalar data
    :param image_name: name of the created image, without extension.
    """
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
            window_size=window_size,
            # title="TPM",
            off_screen=True,
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
            p.update_scalar_bar_range(
                clim=data_layer.color_limits,
                name=data_layer.title,
            )

        # Camera position to zoom to face
        p.camera_position = (0, 3e-4, 1e-3)
        p.camera.zoom(1.1)
        # p.camera.tight(padding=0.05, adjust_render_window=False)
        # print(p.camera)

        output_subdir = Path(output_dir) / f"{name}"
        output_subdir.mkdir(exist_ok=True, parents=True)
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
    for k in track(range(num_steps), description="Create combined images ..."):
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


if __name__ == "__main__":
    # FIXME: make this in an example and add path to resources
    vtk_path_spt = Path("lobule_BCflux.t006.vtk")
    mesh_spt: meshio.Mesh = meshio.read(vtk_path_spt)
    output_path_spt = Path("./raw_spt/")
    output_path_spt.mkdir(exist_ok=True)
    scalars_spt: List[DataLayer] = [
        DataLayer(
            sid="rr_(S)",
            title="Substrate S [mM]",
            colormap="RdBu",
        ),
        DataLayer(
            sid="rr_(P)",
            title="Product P [mM]",
            colormap="RdBu",
        ),
        DataLayer(
            sid="rr_necrosis",
            title="Necrosis",
            colormap="binary",
        ),
    ]
    visualize_data_layers(
        mesh=mesh_spt,
        scalars=scalars_spt,
        output_dir=output_path_spt,
        image_name="example",
    )
