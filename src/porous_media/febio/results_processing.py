"""Process timeseries data in XDMF."""
from pathlib import Path
from typing import Dict

import meshio

from porous_media import DATA_DIR
from porous_media.console import console


def vtks_to_xdmf(vtk_dir: Path, xdmf_path: Path):
    """Convert VTK timesteps to XDMF time course."""
    console.rule(title=f"{xdmf_path}", style="white")

    vtk_paths = sorted(list(vtk_dir.glob("*.vtk")))
    mesh: meshio.Mesh = meshio.read(vtk_paths[0])
    with meshio.xdmf.TimeSeriesWriter(xdmf_path) as writer:
        writer.write_points_cells(mesh.points, mesh.cells)

        for vtk_path in vtk_paths:
            console.print(f"\t{vtk_path}")
            # read timepoint
            with open(vtk_path, "r") as f_vtk:
                f_vtk.readline()  # skip first line
                line = f_vtk.readline()
                t = float(line.strip().split(" ")[-1])

            # read mesh
            mesh = meshio.read(vtk_path)
            writer.write_data(t, point_data=mesh.point_data, cell_data=mesh.cell_data)


# FIXME: get data and interpolate
def interpolate_xdmf(xdmf_in: Path, xdmf_out: Path, times: np.ndarray):
    """Interpolate XDMF."""




# filter vtk by times
times_wanted = np.linspace(0, 600 * 60, 11)  # [s] (21 points in 600 min) # static image
# times_wanted = np.linspace(0, 600*60, 201)  # [s] (21 points in 600 min) # gifs


if __name__ == "__main__":
    # vtk_dirs: Dict[str, Path] = {
    #     "sim_T277": DATA_DIR / "simliva" / "005_T_277_15K_P0__0Pa_t_24h" / "vtk",
    #     "sim_T310": DATA_DIR / "simliva" / "006_T_310_15K_P0__0Pa_t_24h" / "vtk",
    # }
    # for sim_key, vtk_dir in vtk_dirs.items():
    #     vtks_to_xdmf(vtk_dir, xdmf_path=vtk_dir.parent / "results.xdmf")

    interpolate_xdmf(
        xdmf_in=DATA_DIR / "simliva" / "005_T_277_15K_P0__0Pa_t_24h" / "results.xdmf",
        xdmf_out=DATA_DIR / "simliva" / "005_T_277_15K_P0__0Pa_t_24h" / "results_interpolated.xdmf"
    )
