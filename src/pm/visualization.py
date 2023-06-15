"""Visualization of MeshResults."""

from pm.vtk_tools import MeshTimepoint
from pathlib import Path

def visualize_vtk_pyvista(vtk_path: Path) -> None:
    import pyvista as pv
    """Blocking visualization with pyvista."""
    # simply pass the numpy points to the PolyData constructor
    # https://docs.pyvista.org/version/stable/examples/00-load/create-tri-surface.html#sphx-glr-examples-00-load-create-tri-surface-py
    # https://stackoverflow.com/questions/56401123/reading-and-plotting-vtk-file-data-structure-with-python#56482233
    # read the data
    grid = pv.read(vtk_path)

    # plot the data with an automatically created Plotter
    grid.plot(show_scalar_bar=False, show_axes=False)


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

    # block until figure is closed
    mlab.show()


if __name__ == "__main__":
    from pm import EXAMPLE_VTK
    # visualize_vtk_pyvista(vtk_path=EXAMPLE_VTK)
    visualize_vtk_mayavi(vtk_path=EXAMPLE_VTK)