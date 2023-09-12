"""Process timeseries data in XDMF."""
from pathlib import Path
from typing import Dict
import numpy as np
import meshio
import shutil

from porous_media import DATA_DIR
from porous_media.console import console


def vtks_to_xdmf(vtk_dir: Path, xdmf_path: Path) -> None:
    """Convert VTK timesteps to XDMF time course."""
    console.rule(title=f"{xdmf_path}", style="white")

    vtk_paths = sorted(list(vtk_dir.glob("*.vtk")))
    mesh: meshio.Mesh = meshio.read(vtk_paths[0])

    with meshio.xdmf.TimeSeriesWriter(xdmf_path, data_format="HDF") as writer:
        writer.write_points_cells(mesh.points, mesh.cells)

        for vtk_path in vtk_paths[:10]:  # FIXME
            console.print(f"\t{vtk_path}")
            # read timepoint
            with open(vtk_path, "r") as f_vtk:
                f_vtk.readline()  # skip first line
                line = f_vtk.readline()
                t = float(line.strip().split(" ")[-1])

            # read mesh
            mesh = meshio.read(vtk_path)
            writer.write_data(t, point_data=mesh.point_data, cell_data=mesh.cell_data)

    # Fix incorrect *.h5 path
    # https://github.com/nschloe/meshio/pull/1358
    shutil.move(f"{xdmf_path.stem}.h5", str(xdmf_path.parent))


def interpolate_xdmf(xdmf_in: Path, xdmf_out: Path, times_interpolate: np.ndarray):
    """Interpolate XDMF."""
    with meshio.xdmf.TimeSeriesReader(xdmf_in) as reader:
        points, cells = reader.read_points_cells()

        with meshio.xdmf.TimeSeriesWriter(xdmf_out) as writer:
            writer.write_points_cells(points, cells)

            times_data = np.zeros((reader.num_steps,))
            for k in range(reader.num_steps):
                t, _, _ = reader.read_data(k)
                times_data[k] = t

            lower_indices = np.zeros_like(times_interpolate)
            upper_indices = np.zeros_like(times_interpolate)
            for ki, ti in enumerate(times_interpolate):
                lower_indices[ki] = np.argwhere(times_data <= ti)[-1]
                upper_indices[ki] = np.argwhere(times_data >= ti)[0]

            times_lower = [times_data[k] for k in lower_indices]
            times_upper = [times_data[k] for k in upper_indices]
            console.print(f"{times_interpolate}")
            console.print(f"{times_lower}")
            console.print(f"{times_upper}")

            for k, _ in enumerate(times_interpolate):

                _, point_data_lower, cell_data_lower = reader.read_data(lower_indices[k])
                _, point_data_upper, cell_data_upper = reader.read_data(upper_indices[k])

                point_data = {}
                # for key in point_data_lower:
                #     point_data = # FIXME interpolate





if __name__ == "__main__":

    vtk_dirs: Dict[str, Path] = {
        "sim_T277": DATA_DIR / "simliva" / "005_T_277_15K_P0__0Pa_t_24h" / "vtk",
        "sim_T310": DATA_DIR / "simliva" / "006_T_310_15K_P0__0Pa_t_24h" / "vtk",
    }
    for sim_key, vtk_dir in vtk_dirs.items():
        vtks_to_xdmf(vtk_dir, xdmf_path=vtk_dir.parent / "results.xdmf")

    interpolate_xdmf(
        xdmf_in=DATA_DIR / "simliva" / "005_T_277_15K_P0__0Pa_t_24h" / "results.xdmf",
        xdmf_out=DATA_DIR / "simliva" / "005_T_277_15K_P0__0Pa_t_24h" / "results_interpolated1.xdmf",
        times_interpolate=np.linspace(0, 600 * 60, num=11)  # [s] (11 points in 600 min) # static image
    )

    interpolate_xdmf(
        xdmf_in=DATA_DIR / "simliva" / "005_T_277_15K_P0__0Pa_t_24h" / "results.xdmf",
        xdmf_out=DATA_DIR / "simliva" / "005_T_277_15K_P0__0Pa_t_24h" / "results_interpolated2.xdmf",
        times_interpolate=np.linspace(0, 600 * 60, num=250)
    )
