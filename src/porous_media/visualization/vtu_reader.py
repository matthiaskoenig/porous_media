"""Read VTU files for visualization."""

from pathlib import Path
from typing import List

import meshio

from porous_media import RESOURCES_DIR
from porous_media.console import console
from porous_media.data.xdmf_tools import AttributeType
from porous_media.visualization.pyvista_visualization import (
    DataLayer,
    VisualizationSettings,
    visualize_interactive,
)


def vtu_to_xdmf(vtu_path: Path, xdmf_path: Path) -> meshio.Mesh:
    """Convert VTK timesteps to XDMF time course.

    :param overwrite: forces overwrite of files, by default existing files are ignored.
    """
    console.rule(title=f"{xdmf_path}", style="white")

    mesh: meshio.Mesh = meshio.read(vtu_path)

    #
    #
    # if not xdmf_path.parent.exists():
    #     xdmf_path.parent.mkdir(parents=True)
    #
    # console.print(f"{vtk_dir} -> {xdmf_path}")
    # if not overwrite and xdmf_path.exists():
    #     console.print(f"xdmf file exists: {xdmf_path}")
    #     DataLimits.from_xdmf(xdmf_path=xdmf_path, overwrite=overwrite)
    #     return

    return mesh


if __name__ == "__main__":
    vtu_path = RESOURCES_DIR / "vtu" / "1ptestfv-00001.vtu"
    xdmf_path = vtu_path.parent / f"{vtu_path.stem}.xdmf"
    mesh = vtu_to_xdmf(vtu_path=vtu_path, xdmf_path=xdmf_path)

    console.print(mesh)
    console.print(f"{mesh.points}")
    console.print(f"{mesh.cells}")
    console.print(f"{mesh.cell_data}")
    console.print(f"{mesh.point_data}")

    data_layers: List[DataLayer] = [
        DataLayer(
            sid="p",
            title="pressure p",
            colormap="magma",
            viz_type=AttributeType.SCALAR,
        ),
        DataLayer(
            sid="velocity_liq (m/s)",
            title="velocity_liq (m/s)",
            colormap="magma",
            viz_type=AttributeType.VECTOR,
        ),
        DataLayer(
            sid="process rank",
            title="process rank",
            colormap="magma",
            viz_type=AttributeType.SCALAR,
        ),
    ]

    visualize_interactive(
        mesh, data_layer=data_layers[2], visualization_settings=VisualizationSettings()
    )
