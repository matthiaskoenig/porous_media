"""Mayavi application."""
from pathlib import Path

def visualize_vtk_mayavi(vtk_path: Path):

    from mayavi import mlab
    from mayavi.modules.surface import Surface

    # visualize mesh

    # create a new figure, grab the engine that's created with it
    fig = mlab.figure()
    engine = mlab.get_engine()

    # open the vtk file, let mayavi figure it all out
    vtk_file_reader = engine.open(str(vtk_path))

    # plot surface corresponding to the data
    surface = Surface()
    engine.add_filter(surface, vtk_file_reader)

    scene = engine.scenes[0]
    scene.scene.background = (1.0, 1.0, 1.0)
    scene.scene.z_minus_view()
    scene.scene.show_axes = False
    module_manager = engine.scenes[0].children[0].children[0]
    module_manager.scalar_lut_manager.lut_mode = 'RdBu'

    # block until figure is closed
    mlab.show()


if __name__ == "__main__":
    from pm import EXAMPLE_VTK
    # visualize_vtk_pyvista(vtk_path=EXAMPLE_VTK)
    visualize_vtk_mayavi(vtk_path=EXAMPLE_VTK)
