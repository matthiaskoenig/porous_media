
from pm.console import console
from pm.vtk_tools import MeshTimepoint, read_mesh_timepoint_from_vtk
from pm.visualization import visualize_mesh_pyvista


if __name__ == "__main__":
    from pm import EXAMPLE_VTK
    mesh_results: MeshTimepoint = read_mesh_timepoint_from_vtk(EXAMPLE_VTK)

    # visualize_mesh_pyvista(mesh_results.mesh)




