from pathlib import Path

import meshio
import pyvista as pv

if __name__ == "__main__":
    vtk_example = Path(__file__).parent.parent.parent / "data" / "sim001_perf1_BC1.t025.vtk"
    mesh = meshio.read(vtk_example)
    print(mesh)
    print("points:", mesh.points)
    print("cells:", mesh.cells)

    print("data:", mesh.point_data.keys())
    print("data:", mesh.cell_data.keys())
    # mesh.points, mesh.cells, ...

    # simply pass the numpy points to the PolyData constructor
    # https://docs.pyvista.org/version/stable/examples/00-load/create-tri-surface.html#sphx-glr-examples-00-load-create-tri-surface-py
    cloud = pv.PolyData(mesh.points)
    cloud.plot(point_size=15)
    surf = cloud.delaunay_2d()
    surf.plot(show_edges=True)
