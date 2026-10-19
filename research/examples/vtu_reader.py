"""Archived example; configure research paths before running."""

from porous_media import RESOURCES_DIR
from porous_media.console import console
from porous_media.data.xdmf_tools import AttributeType
from porous_media.visualization.pyvista_visualization import (
    DataLayer,
    VisualizationSettings,
    visualize_interactive,
)
from porous_media.visualization.vtu_reader import vtu_to_xdmf
from research import RESULTS_DIR

if __name__ == "__main__":
    vtu_path = RESOURCES_DIR / "vtu" / "1ptestfv-00001.vtu"
    xdmf_path = RESULTS_DIR / "vtu" / f"{vtu_path.stem}.xdmf"
    mesh = vtu_to_xdmf(vtu_path=vtu_path, xdmf_path=xdmf_path)

    console.print(mesh)
    console.print(f"{mesh.points}")
    console.print(f"{mesh.cells}")
    console.print(f"{mesh.cell_data}")
    console.print(f"{mesh.point_data}")

    data_layers: list[DataLayer] = [
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
