"""Visualization with pyvista."""
import tempfile
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import meshio
import numpy as np
import pyvista as pv
from rich.progress import track

from porous_media.console import console
from porous_media.mesh.mesh_tools import MeshTimepoint, mesh_to_xdmf
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
class Scalar:
    """Scalar information for plotting."""

    sid: str
    title: str
    colormap: str
    color_limits: Optional[List[float]] = None
    scalar_bar: bool = True


def calculate_value_ranges(
    xdmf_path: Path, scalars: List[Scalar]
) -> Dict[str, List[float]]:
    """Calculate the limits/ranges for given scalars and adds to the scalars.

    Iterates over all timepoints of a given scalar to figure out the ranges of the data.
    """
    if not xdmf_path:
        raise IOError(f"xdmf_path does not exist: {xdmf_path}")

    limits: Dict[str, List[float]] = {}
    with meshio.xdmf.TimeSeriesReader(xdmf_path) as reader:
        points, cells = reader.read_points_cells()
        # iterate over time points
        for k in range(reader.num_steps):
            t, point_data, cell_data = reader.read_data(k)

            for scalar in scalars:
                scalar_limits: List[float] = [np.NaN, np.NaN]

                # process cell data
                sid = scalar.sid
                try:
                    data = cell_data[scalar.sid]
                except KeyError as err:
                    console.print(f"Exising keys: {cell_data.keys()}")
                    raise err

                data = np.dstack(data)
                data = data.squeeze()

                # new min, max
                dmin = data.min()
                dmax = data.max()
                if not scalar_limits:
                    scalar_limits = [dmin, dmax]
                else:
                    if scalar_limits[0] > dmin:
                        scalar_limits[0] = dmin
                    if scalar_limits[1] < dmax:
                        scalar_limits[1] = dmax

                scalar.color_limits = scalar_limits
                limits[sid] = scalar_limits

    return limits


def visualize_scalars_timecourse(
    xdmf_path: Path,
    output_dir: Path,
    scalars: List[Scalar],
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
        for k in track(range(tnum), description="Create panels for scalars..."):
            t, point_data, cell_data = reader.read_data(k)

            # Create mesh with single data point
            mesh: meshio = meshio.Mesh(
                points=points, cells=cells, cell_data=cell_data, point_data=point_data
            )
            visualize_scalars(
                mesh=mesh,
                scalars=scalars,
                image_name=f"sim_{k:05d}",
                output_dir=output_dir / "panels",
                window_size=window_size,
            )


def visualize_scalars(
    mesh: meshio.Mesh,
    scalars: List[Scalar],
    output_dir: Path,
    image_name: str,
    window_size: Tuple[float, float] = (1000, 1000),
) -> None:
    """Visualize geometry with pyvista.

    :param mesh: mesh with single time point scalar data
    :param image_name: name of the created image, without extension.
    """

    # create VTK from mesh to read for pyvista
    # FIXME: better handling of grid
    xdmf_tmp = tempfile.NamedTemporaryFile(suffix=".xdmf")
    mesh_to_xdmf(m=mesh, xdmf_path=Path(xdmf_tmp.name), test_read=False)
    grid = pv.read(xdmf_tmp.name)

    # deactivate active sets
    grid.set_active_tensors(None)
    grid.set_active_scalars(None)
    grid.set_active_vectors(None)

    # visualize scalars
    scalars_dict: Dict[str, Scalar] = {s.sid: s for s in scalars}
    for scalar_id, scalar in scalars_dict.items():
        p = pv.Plotter(
            window_size=window_size,
            title="TPM lobulus",
            off_screen=True,
        )

        scalar = scalars_dict[scalar_id]
        grid.set_active_scalars(name=scalar_id)

        actor = p.add_mesh(
            grid,
            show_edges=True,
            render_points_as_spheres=True,
            point_size=3,
            show_vertices=True,
            line_width=1.0,
            cmap=scalar.colormap,
            show_scalar_bar=False,
            edge_color="darkgray",
            # specular=0.5, specular_power=15,
            clim=scalar.color_limits,
        )
        if scalar.scalar_bar:
            p.add_scalar_bar(
                title=scalar.title,
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
                clim=scalar.color_limits,
                name=scalar.title,
            )

        # Camera position to zoom to face
        p.camera_position = (0, 3e-4, 1e-3)
        p.camera.zoom(1.1)
        # p.camera.tight(padding=0.05, adjust_render_window=False)
        # print(p.camera)

        output_subdir = Path(output_dir) / f"{scalar_id}"
        output_subdir.mkdir(exist_ok=True, parents=True)
        p.show(
            # cpos=top_view,
            screenshot=Path(output_subdir / f"{image_name}.png")
        )


def create_combined_images(
    xdmf_path: Path,
    output_dir: Path,
    scalars_selection: List[str],
    direction: str,
    ncols: Optional[int] = None,
    nrows: Optional[int] = None,
) -> List[Path]:
    """Create combined images for all time points."""
    image_dir: Path = output_dir / direction
    image_dir.mkdir(parents=True, exist_ok=True)

    all_images: List[Path] = []
    with meshio.xdmf.TimeSeriesReader(xdmf_path) as reader:
        _, _ = reader.read_points_cells()
        for k in track(
            range(reader.num_steps), description="Create combined images ..."
        ):
            images: List[Path] = []
            for scalar_id in scalars_selection:
                img_path = output_dir / "panels" / scalar_id / f"sim_{k:05d}.png"
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
    scalars_spt: List[Scalar] = [
        Scalar(
            sid="rr_(S)",
            title="Substrate S [mM]",
            colormap="RdBu",
        ),
        Scalar(
            sid="rr_(P)",
            title="Product P [mM]",
            colormap="RdBu",
        ),
        Scalar(
            sid="rr_necrosis",
            title="Necrosis",
            colormap="binary",
        ),
    ]
    visualize_scalars(
        mesh=mesh_spt,
        scalars=scalars_spt,
        output_dir=output_path_spt,
        image_name="example",
    )
