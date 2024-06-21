"""Add position variable and zonation patterns to mesh."""

from porous_media.mesh.mesh_tools import mesh_to_vtk, mesh_to_xdmf
from porous_media.mesh.mesh_zonation import ZonatedMesh
import meshio

# FIXME: require: "element_volume_point_TPM"

from pathlib import Path

results_path: Path = Path(__file__).parent
vtk_path = results_path / "mesh_no_pressure.vtk"

# create zonated mesh
zm = ZonatedMesh()
m: meshio.Mesh = zm.create_zonated_mesh_from_vtk(vtk_path=vtk_path)

# serialize mesh with results
xdmf_protein_path = vtk_path.parent / f"{vtk_path.stem}_protein.xdmf"
mesh_to_xdmf(m=m, xdmf_path=xdmf_protein_path)
vtk_protein_path = results_path / f"{vtk_path.stem}_protein.vtk"
mesh_to_vtk(m=m, vtk_path=vtk_protein_path)
