"""Visualization with pyvista."""
import tempfile
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import meshio
import numpy as np
import pyvista as pv

from porous_media.console import console
from porous_media.mesh.mesh_tools import mesh_to_xdmf, MeshTimepoint
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
    cmap: str
    color_limits: Optional[Tuple[float, float]] = None


def calculate_value_ranges(xdmf_path: Path, scalars: List[Scalar]) -> Dict[str, List[float]]:
    """Calculate the limits/ranges for given scalars and adds to the scalars."""

    limits = {}
    with meshio.xdmf.TimeSeriesReader(xdmf_path) as reader:
        points, cells = reader.read_points_cells()

        # iterate over time points
        for k in range(reader.num_steps):
            t, point_data, cell_data = reader.read_data(k)

            for scalar in scalars:
                scalar_limits = None

                # process cell data
                sid = scalar.sid
                data = cell_data[scalar.sid]
                data = np.dstack(data)
                data = data.squeeze()
                # console.print(f"{data=}")
                # console.print(f"{type(data)}")

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

                scalar.color_limits=scalar_limits
                limits[sid] = scalar_limits

    console.print(limits)
    return limits

def visualize_scalars(
    mesh: meshio.Mesh,
    results_path: Path,
    image_name: str,
    drange_type: DataRangeType = DataRangeType.LOCAL,
) -> None:
    """Visualize zonation patterns.

    :param image_name: Name for the generated image in the output path.
    """

    # calculate the data ranges
    dmins: Dict[str, float] = {}
    dmaxs: Dict[str, float] = {}
    for key in mesh.cell_data:
        data = mesh.cell_data[key][0]
        # min, max
        dmins[key] = data.min()
        dmaxs[key] = data.max()

    # handle the case of global data ranges
    if drange_type == DataRangeType.GLOBAL:
        dmin_global = min(dmins.values())
        dmax_global = max(dmaxs.values())
        for key in mesh.cell_data:
            dmins[key] = dmin_global
            dmaxs[key] = dmax_global

    scalars = {}
    for key in mesh.cell_data:
        if key.startswith("pattern__") or key == "cell_type":
            pattern = key.split("__")[-1]
            scalars[key] = {
                "title": f"{pattern.upper()} [-]",
                "cmap": "RdBu",  # "cmap": "Blues", # FIXME: better colormap
                "clim": (dmins[key], dmaxs[key]),
            }

    # create raw images
    visualize_lobulus_vtk(
        mesh=mesh, scalars=scalars, output_dir=results_path, image_name=image_name
    )
    # combine images
    images: List[Path] = []
    for scalar in scalars:
        img_path = results_path / scalar / f"{image_name}.png"
        images.append(img_path)

    image: Path = results_path / f"{image_name}.png"
    merge_images(paths=images, direction="horizontal", output_path=image)
    console.print(f"Image created: file://{image}")


def visualize_lobulus_vtk(
    mesh: meshio.Mesh,
    scalars: Dict,
    output_dir: Path,
    image_name: str,
    window_size: Tuple[float, float] = (1000, 1000),
    scalar_bar: bool = True,
) -> None:
    """Visualize single lobulus time point with pyvista.

    :param image_name: name of the created image, without extension.
    """

    # create VTK from mesh to read for pyvista
    xdmf_tmp = tempfile.NamedTemporaryFile(suffix=".xdmf")
    mesh_to_xdmf(m=mesh, xdmf_path=Path(xdmf_tmp.name), test_read=False)
    grid = pv.read(xdmf_tmp.name)

    # deactivate active sets
    # grid.set_active_tensors(None)
    grid.set_active_scalars(None)
    # grid.set_active_vectors(None)

    # plot the data with an automatically created Plotter
    # grid.plot(show_scalar_bar=False, show_axes=False)

    for scalar_id in scalars:
        p = pv.Plotter(
            window_size=window_size,
            title="TPM lobulus",
            # shape=(1, n_scalars), border=False,
            off_screen=True,
        )

        scalar_info = scalars[scalar_id]
        grid.set_active_scalars(name=scalar_id)

        actor = p.add_mesh(
            grid,
            show_edges=True,
            render_points_as_spheres=True,
            point_size=3,
            show_vertices=True,
            line_width=1.0,
            cmap=scalar_info["cmap"],
            show_scalar_bar=False,
            edge_color="darkgray",
            # specular=0.5, specular_power=15,
            clim=scalar_info["clim"],
        )
        if scalar_bar:
            p.add_scalar_bar(
                title=scalar_info["title"],
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
                clim=scalar_info["clim"], name=scalar_info["title"]
            )

        # Camera position to zoom to face
        p.camera_position = (0, 3e-4, 1e-3)
        p.camera.zoom(1.1)
        # p.camera.tight(padding=0.05, adjust_render_window=False)
        # print(p.camera)

        output_subdir = Path(output_dir) / f"{scalar_id}"
        output_subdir.mkdir(exist_ok=True)
        p.show(
            # cpos=top_view,
            screenshot=Path(output_subdir / f"{image_name}.png")
        )


def visualize_panels(
    xdmf_path: Path, output_path: Path, scalars: List[Scalar]
) -> None:
    """Visualize the panels."""


    # FIXME: update using the xdmf
    # Create figures for all simulation timepoints of relevance (skip some results)
    for sim_key, vtk_paths in vtks.items():
        panels_path = output_path / sim_key / "panels"
        panels_path.mkdir(exist_ok=True, parents=True)

        # Create all figures for all variables
        for vtk_path in vtk_paths:
            show_scalar_bar = True
            # if k == 0 or k == (len(vtk_paths)-1):
            #     show_scalar_bar = True
            mesh: meshio.Mesh = meshio.read(vtk_path)
            console.print(mesh)
            visualize_lobulus_vtk(
                mesh=mesh,
                scalars=scalars,
                image_name=vtk_path.stem,
                output_dir=panels_path,
                window_size=(600, 600),
                scalar_bar=show_scalar_bar,
            )

if __name__ == "__main__":
    # TODO: get the clims for the scalars from all vtks or global settings

    # FIXME: make this in an example and add path to resources
    vtk_path_spt = Path("lobule_BCflux.t006.vtk")
    mesh_spt: meshio.Mesh = meshio.read(vtk_path_spt)
    output_path_spt = Path("./raw_spt/")
    output_path_spt.mkdir(exist_ok=True)
    scalars_spt = {
        "rr_(S)": {
            "title": "Substrate S [mM]",
            "cmap": "RdBu",
            # "clim": (0.529, 0.530),
        },
        "rr_(P)": {
            "title": "Product P [mM]",
            "cmap": "RdBu",
            # "clim": (0.00829, 0.00897),
        },
        "rr_necrosis": {
            "title": "Necrosis",
            "cmap": "binary",
            # "clim": (0.0, 1.0),
        },
    }
    visualize_lobulus_vtk(
        mesh=mesh_spt,
        scalars=scalars_spt,
        output_dir=output_path_spt,
        image_name="example",
    )
