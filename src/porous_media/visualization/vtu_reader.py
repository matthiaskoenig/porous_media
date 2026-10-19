"""Read VTU files for visualization."""

from pathlib import Path

import meshio

from porous_media.console import console
from porous_media.mesh.mesh_tools import mesh_to_xdmf


def vtu_to_xdmf(
    vtu_path: Path, xdmf_path: Path, *, overwrite: bool = False
) -> meshio.Mesh:
    """Convert a static VTU mesh to XDMF/HDF and return the mesh.

    Args:
        vtu_path: Input VTU mesh path.
        xdmf_path: Destination XDMF path. Its HDF sidecar stays beside it.
        overwrite: Explicitly allow replacement of existing regular output files.

    Returns:
        The mesh read from the VTU file.
    """
    console.rule(title=f"{xdmf_path}", style="white")
    if vtu_path.resolve() == xdmf_path.resolve():
        raise ValueError("Input and output paths must differ")
    for output in (xdmf_path, xdmf_path.with_suffix(".h5")):
        if output.is_symlink():
            raise ValueError("XDMF outputs must not follow symbolic links")
        if output.exists() and not overwrite:
            raise FileExistsError(f"Output already exists: {output}")

    mesh: meshio.Mesh = meshio.read(vtu_path)
    xdmf_path.parent.mkdir(parents=True, exist_ok=True)
    mesh_to_xdmf(mesh, xdmf_path)
    return mesh
