"""Mayavi script.

**datasets**:
https://docs.enthought.com/mayavi/mayavi/data.html

All VTK arrays, whether it be for data or position, are exposed as (n, 3) numpy arrays
for 3D components, and flat (n, ) array for 1D components. The index vary in the
opposite order as numpy: z first, y and then x. Thus to go from a 3D numpy array to
the corresponding flatten VTK array, the operation is:

vtk_array = numpy_array.T.ravel()

**filters**
transform data, e.g. PointToCellData or CellToPointData

SelectOutput !!! A filter that allows a user to select one among several of the
outputs of a given input. This is typically very useful for a multi-block data source.

"""
from pathlib import Path

import numpy as np
from mayavi import mlab
from pm.mesh_tools import MeshTimepoint


def visualize_lobulus(mtp: MeshTimepoint):
    # create figure
    mlab.figure("Fig1", bgcolor=(1.0, 1.0, 1.0), fgcolor=(0.0, 0.0, 0.0), size=(800, 800))


    # # Example data
    # from numpy import pi, sin, cos, mgrid
    # dphi, dtheta = pi / 250.0, pi / 250.0
    # [phi, theta] = mgrid[0:pi + dphi * 1.5:dphi, 0:2 * pi + dtheta * 1.5:dtheta]
    # m0 = 4;
    # m1 = 3;
    # m2 = 2;
    # m3 = 3;
    # m4 = 6;
    # m5 = 2;
    # m6 = 6;
    # m7 = 4;
    # r = sin(m0 * phi) ** m1 + cos(m2 * phi) ** m3 + sin(m4 * theta) ** m5 + cos(m6 * theta) ** m7
    # x = r * sin(phi) * cos(theta)
    # y = r * cos(phi)
    # z = r * sin(phi) * sin(theta)
    #
    # # render data with colormap
    # m = mlab.mesh(
    #     x, y, z,
    #     colormap='RdBu',
    # )

    # To retrieve the VTK datasets feeding in an arbitrary object, the mlab
    # function pipeline.get_vtk_src() may be useful.
    # mlab.pipeline.get_vtk_src(m)


    # TODO: load meshio data from VTK into the mesh
    # x, y, z are arrays giving the positions of the vertices of the surface.
    # triangles is a list of triplets (or an array) list the vertices in each triangle.
    mesh = mtp.mesh
    console.print("points:", mesh.points)
    npoints = len(mesh.points)
    x = np.zeros(npoints)
    y = np.zeros(npoints)
    z = np.zeros(npoints)

    for k, point in enumerate(mesh.points):
        x[k] = point[0]
        y[k] = point[1]
        z[k] = point[2]

    data = x + y
    print(x.shape)

    src = mlab.pipeline.scalar_scatter(
        x, y, z,
        data
    )
    pts = mlab.pipeline.glyph(
        src,
        scale_mode='none', scale_factor=.0001,
        color=(0, 0, 0)
    )

    console.print(z)
    # reduce to a 2d plane via the z coordinate

    # structured grid
    # from tvtk.api import tvtk
    # sgrid = tvtk.StructuredGrid(dimensions=(dims[1], dims[0], dims[2]))
    # sgrid.points = pts
    # s = random.random((dims[0] * dims[1] * dims[2]))
    # sgrid.point_data.scalars = ravel(s.copy())
    # sgrid.point_data.scalars.name = 'scalars'

    # unstructured grid


    console.print("cells:", mesh.cells[0])
    from meshio import CellBlock
    cell_block: CellBlock = mesh.cells[0]
    console.print(cell_block.data)

    # get the triangles from the wedge

    # get triangle from wedge:


    # triangles = [] figure out the triangles

    # 'rr_(S)': (563,),  # cell data
    d = mesh.cell_data['rr_(S)']

    # mlab.triangular_mesh(
    #     x, y, z,
    #     # triangles
    #     # scalars=d,
    #     representation="surface",
    #     colormap='RdBu',
    # )
    # representation: the representation type used for the surface. Must be ‘surface’ or ‘wireframe’ or ‘points’ or ‘mesh’ or ‘fancymesh’. Default: surface


    # add colorbar
    mlab.colorbar(title="z data", orientation="vertical", label_fmt="%.1f", nb_labels=5)

    # set camera position
    # TODO: set top level view
    mlab.view(
        azimuth=None,
        elevation=None,
        distance=None,
        focalpoint=None,
        roll=None, reset_roll=True
    )

    # screenshot
    mlab.savefig("test.png", magnification=2)

    mlab.show()




if __name__ == "__main__":

    from pm.console import console


    for vtk_path in [
        "lobule_BCflux.t006.vtk",
        "t2793.vtk",
    ]:
        # TODO: load mesh data with content
        mtp: MeshTimepoint = MeshTimepoint.from_vtk(vtk_path=vtk_path, show=True)
        console.print("cell_type")
        console.print(mtp.cell_data["cell_type"])
        console.print(mtp)
        console.rule()

        # # we want to read VTK in an UnstructuredGrid object
        # from mayavi.sources import vtk_file_reader
        # r = vtk_file_reader.VTKFileReader()
        # r.initialize(str(vtk_path))
        #
        # r.
        #
        # # use the source
        #
        #
        # # open the vtk file, let mayavi figure it all out
        # vtk_file_reader = engine.open(file_name)
        #
        #
        # visualize_lobulus(mtp=mtp)

        # TODO: get the x, y, z data for the mesh