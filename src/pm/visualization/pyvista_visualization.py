"""Visualization with pyvista"""
from enum import Enum
from pathlib import Path
from typing import Dict, List

import meshio
import pyvista as pv

from pm.console import console
from pm.mesh.mesh_tools import mesh_to_vtk
from pm.visualization.image_manipulation import merge_images


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
    window_size=(1000, 1000),
    scalar_bar: bool = True,
):
    """Visualize single lobulus time point with pyvista.

    :param image_name: name of the created image, without extension.
    """

    # create VTK from mesh to read for pyvista
    import tempfile

    vtk_tmp = tempfile.NamedTemporaryFile(suffix=".vtk")
    mesh_to_vtk(mesh, vtk_tmp.name, test_read=False)

    # read the data
    grid = pv.read(vtk_tmp.name)
    # console.print(grid)
    # console.print(grid.cell_data)

    # deactivate active sets
    # grid.set_active_tensors(None)
    grid.set_active_scalars(None)
    # grid.set_active_vectors(None)

    # plot the data with an automatically created Plotter
    # grid.plot(show_scalar_bar=False, show_axes=False)

    n_scalars = len(scalars)
    for k, scalar_id in enumerate(scalars):
        p = pv.Plotter(
            window_size=window_size,
            title="TPM lobulus",
            # shape=(1, n_scalars), border=False,
            off_screen=True,
        )

        scalar_info = scalars[scalar_id]

        # p.subplot(0, k)
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


if __name__ == "__main__":
    # TODO: get the clims for the scalars from all vtks or global settings

    # TODO: process all files of simulation

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
        mesh=mesh_spt, scalars=scalars_spt, output_dir=output_path_spt
    )

    console.rule(style="white")

    vtk_path_iri = Path("t2793.vtk")
    output_path_iri = Path("./raw_iri/")
    output_path_iri.mkdir(exist_ok=True)
    scalars_iri = {
        "ATP": {
            "title": "ATP [mM]",
            "cmap": "RdBu",
            # "clim": (0.529, 0.530),
        },
        "ADP": {
            "title": "ADP [mM]",
            "cmap": "RdBu",
            # "clim": (0.00829, 0.00897),
        },
        "necrosis": {
            "title": "Necrosis",
            "cmap": "binary",
            # "clim": (0.0, 1.0),
        },
    }
    visualize_lobulus_vtk(
        vtk_path=vtk_path_iri, scalars=scalars_iri, output_dir=output_path_iri
    )
