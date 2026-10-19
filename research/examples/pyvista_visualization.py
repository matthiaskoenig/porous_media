"""Archived example; configure research paths before running."""

from porous_media import RESOURCES_DIR
from porous_media.console import console
from porous_media.data.xdmf_tools import AttributeType, XDMFInfo, vtks_to_xdmf
from porous_media.visualization.pyvista_visualization import (
    DataLayer,
    VisualizationSettings,
    visualize_interactive,
    xdmf_to_mesh,
)
from research import RESULTS_DIR

if __name__ == "__main__":
    # FIXME: make this in an example and add path to resources

    # Simple test that all variables are read

    # vtk_dir = RESOURCES_DIR / "vtk" / "vtk_single"
    vtk_dir = RESOURCES_DIR / "vtk" / "vtk_timecourse"
    xdmf_path = RESULTS_DIR / "vtk_test.xdmf"

    vtks_to_xdmf(vtk_dir, xdmf_path=xdmf_path, overwrite=True)
    xdmf_info: XDMFInfo = XDMFInfo.from_path(xdmf_path)
    console.print(xdmf_info)

    data_layers: list[DataLayer] = [
        DataLayer(
            sid="displacement",
            title="displacement",
            colormap="magma",
            viz_type=AttributeType.VECTOR,
        ),
        DataLayer(
            sid="fluid_flux_TPM",
            title="fluid flux [m/s]",
            colormap="magma",
            viz_type=AttributeType.VECTOR,
        ),
        DataLayer(
            sid="pressure",
            title="Pressure",
            colormap="magma",
            viz_type=AttributeType.SCALAR,
        ),
        DataLayer(
            sid="effective_fluid_pressure_TPM",
            title="Effective fluid pressure [Pa]",
            colormap="magma",
            viz_type=AttributeType.SCALAR,
        ),
        DataLayer(
            sid="rr_(S)",
            title="Substrate S [mM]",
            colormap="magma",
            viz_type=AttributeType.SCALAR,
        ),
    ]
    # visualize_datalayers_timecourse(
    #     xdmf_path=xdmf_path,
    #     data_layers=data_layers,
    #     output_dir=RESULTS_DIR / "vtk_test_images"
    # )

    mesh = xdmf_to_mesh(xdmf_path, k=1)
    visualize_interactive(
        mesh, data_layer=data_layers[2], visualization_settings=VisualizationSettings()
    )
    # FIXME: support calculation of new variables for visualization
