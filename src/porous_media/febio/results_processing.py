"""Process timeseries data in XDMF."""
import os
from pathlib import Path
from typing import Dict, List
from rich.progress import track
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

        for k in track(range(len(vtk_paths)), description="Processing VTKs ..."):
            vtk_path = vtk_paths[k]
            # console.print(f"\t{vtk_path}")
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
    h5_path = xdmf_path.parent / f"{xdmf_path.stem}.h5"
    if h5_path.exists():
        os.remove(h5_path)
    shutil.move(f"{xdmf_path.stem}.h5", str(xdmf_path.parent))


def interpolate_xdmf(xdmf_in: Path, xdmf_out: Path, times_interpolate: np.ndarray):
    """Interpolate XDMF."""
    console.rule(title=f"Interpolate {xdmf_in}", style="white")
    with meshio.xdmf.TimeSeriesReader(xdmf_in) as reader:
        points, cells = reader.read_points_cells()

        with meshio.xdmf.TimeSeriesWriter(xdmf_out) as writer:
            writer.write_points_cells(points, cells)

            tnum = reader.num_steps
            times_data = np.zeros((tnum,))
            for k in track(range(tnum), description="Calculate timepoints ..."):
                t, _, _ = reader.read_data(k)
                times_data[k] = t

            if times_interpolate[0] < times_data[0]:
                raise ValueError(f"Lower interpolation range outside of data: {times_interpolate[0]} < {times_data[0]}")
            if times_interpolate[-1] > times_data[-1]:
                raise ValueError(f"Upper interpolation range outside of data: {times_interpolate[-1]} > {times_data[-1]}")

            lower_indices = np.zeros_like(times_interpolate, dtype=int)
            upper_indices = np.zeros_like(times_interpolate, dtype=int)
            for ki, ti in enumerate(times_interpolate):
                lower_indices[ki] = np.argwhere(times_data <= ti).flatten()[-1]
                upper_indices[ki] = np.argwhere(times_data >= ti).flatten()[0]

            # interpolate data for all data points
            for k in track(range(len(times_interpolate)), description="Interpolate data ..."):

                t_interpolate = times_interpolate[k]
                idx_low = lower_indices[k]
                idx_up = upper_indices[k]

                # interpolate all numpy matrices
                t_low, point_data_low, cell_data_low = reader.read_data(idx_low)
                t_up, point_data_up, cell_data_up = reader.read_data(idx_up)
                if np.isclose(t_up, t_low):
                    f_low = 0.0
                    f_up = 1.0
                else:
                    f_low: float = (t_interpolate - t_low) / (t_up - t_low)
                    f_up: float = (t_up - t_interpolate) / (t_up - t_low)

                # interpolate point data
                point_data = {}
                for key in point_data_low:
                    point_data[key] = (f_low * point_data_low[key]) + (f_up * point_data_up[key])

                # interpolate cell data
                cell_data = {}
                for key in cell_data_low:
                    # console.print(f"{cell_data_low[key]}")
                    # process the list
                    cell_data[key] = [(f_low * cell_data_low[key][k]) + (f_up * cell_data_low[key][k]) for k in range(len(cell_data_low[key]))]

                # write interpolation point
                writer.write_data(t_interpolate, point_data=point_data, cell_data=cell_data)

        # Fix incorrect *.h5 path
        # https://github.com/nschloe/meshio/pull/1358
        h5_path = xdmf_out.parent / f"{xdmf_out.stem}.h5"
        if h5_path.exists():
            os.remove(h5_path)
        shutil.move(f"{xdmf_out.stem}.h5", str(xdmf_out.parent))

        console.print(f"Interpolated data: {xdmf_out}")


if __name__ == "__main__":
    sim_flux_path: Path = DATA_DIR / "simliva" / "iri_flux_study_0" / "vtk"
    sim_277k_path: Path = DATA_DIR / "simliva" / "005_T_277_15K_P0__0Pa_t_24h" / "vtk"
    sim_310k_path: Path = DATA_DIR / "simliva" / "006_T_310_15K_P0__0Pa_t_24h" / "vtk"

    # interpolate_xdmf(
    #     xdmf_in=DATA_DIR / "simliva" / "005_T_277_15K_P0__0Pa_t_24h" / "results.xdmf",
    #     xdmf_out=DATA_DIR / "simliva" / "005_T_277_15K_P0__0Pa_t_24h" / "results_interpolated_11.xdmf",
    #     times_interpolate=np.linspace(0, 600 * 60, num=11)  # [s] (11 points in 600 min) # static image
    # )
    #
    # interpolate_xdmf(
    #     xdmf_in=DATA_DIR / "simliva" / "005_T_277_15K_P0__0Pa_t_24h" / "results.xdmf",
    #     xdmf_out=DATA_DIR / "simliva" / "005_T_277_15K_P0__0Pa_t_24h" / "results_interpolated_11.xdmf",
    #     times_interpolate=np.linspace(0, 600 * 60, num=11)
    # )
    #
    # interpolate_xdmf(
    #     xdmf_in=sim_flux_path.parent / "results.xdmf",
    #     xdmf_out=sim_flux_path.parent / "results_interpolated.xdmf",
    #     times_interpolate=np.linspace(0, 10000, num=200)
    # )

    nums = [10, 200]
    for num in nums:
        interpolate_xdmf(
            xdmf_in=sim_flux_path.parent / "results.xdmf",
            xdmf_out=sim_flux_path.parent / f"results_interpolated_{num}.xdmf",
            times_interpolate=np.linspace(0, 10000, num=num)
        )

