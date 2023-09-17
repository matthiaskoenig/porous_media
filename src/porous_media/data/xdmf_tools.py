"""Process timeseries data in XDMF.

https://www.xdmf.org/index.php/XDMF_Model_and_Format

"""
from __future__ import annotations

import os
import shutil
from copy import deepcopy
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Dict, List, Any, Tuple

import meshio
import numpy as np
from rich.progress import track

from meshio.xdmf.common import attribute_type
from meshio._common import raw_from_cell_data
from porous_media.console import console
from porous_media.log import get_logger

logger = get_logger(__name__)


@dataclass
class XDMFInformation:
    """Information about a given timecourse xdmf.

    This includes information about the timepoints, scalars, vectors, ... .
    This is often required as base information before processing things.

    - figure out time points and range
    - figure out data limits
    - figure out scalars, vectors, ...

    Preprocessing of all the important information for the visualization later on.
    """

    path: Path
    tstart: float
    tend: float
    num_steps: int
    points: List[Any]
    cells: List[Any]
    point_data: Dict[str, str]
    cell_data: Dict[str, str]

    @staticmethod
    def from_path(xdmf_path: Path) -> XDMFInformation:
        """Create XDMFInformation from xdmf path."""
        if not xdmf_path:
            raise IOError(f"xdmf_path does not exist: {xdmf_path}")

        with meshio.xdmf.TimeSeriesReader(xdmf_path) as reader:
            points, cells = reader.read_points_cells()

            # process most information based on initial data and points and cells
            tstart, point_data, cell_data = reader.read_data(0)

            point_data_info = {}
            for name, data in point_data.items():
                point_data_info[name] = attribute_type(data)

            cell_data_info = {}
            raw = raw_from_cell_data(cell_data)
            for name, data in raw.items():
                cell_data_info[name] = attribute_type(data)

            # get end time
            tend, _, _ = reader.read_data(reader.num_steps-1)

        xdmf_info = XDMFInformation(
            path=xdmf_path,
            tstart=tstart,
            tend=tend,
            num_steps=reader.num_steps,
            points=points,
            cells=cells,
            point_data={},
            cell_data=cell_data_info,
        )
        return xdmf_info

@dataclass
class DataLimits:
    """Limits of numerical data."""

    limits: Dict[str, Tuple[float, float]]

    @staticmethod
    def merge_limits(data_limits_list: List[DataLimits]) -> DataLimits:
        """Merge the limits from multiple data limits."""


        limits = deepcopy(data_limits_list[0].limits)
        for k, dlim in enumerate(data_limits_list):
            if k==0:
                continue
            else:
                for name, lims_k in enumerate(dlim.limits):
                    lims = limits[name]
                    if lims_k[0] < lims[0]:
                        FIXME




    @staticmethod
    def from_xdmf(xdmf_path: Path) -> DataLimits:
        """Helper for calculating all data limits.

        This takes some time because it requires iteration over the complete dataset.
        Should support subsets and operations such as merging of limits from multiple
        simulations.
        """
        # read information
        xdmf_info = XDMFInformation.from_path(xdmf_path)

        point_limits: Dict[str, List[float]] = {}
        cell_limits: Dict[str, List[float]] = {}
        with meshio.xdmf.TimeSeriesReader(xdmf_path) as reader:
            _, _ = reader.read_points_cells()

            for k in track(range(reader.num_steps), description="Calculate data limits ..."):
                t, point_data, cell_data = reader.read_data(k)

                # point data limits
                for name, attribute_type in xdmf_info.point_data.items():
                    # min and max from data
                    data = point_data[name]
                    data = np.dstack(data)
                    data = data.squeeze()
                    dmin = data.min()
                    dmax = data.max()

                    limits: Tuple[float, float]
                    if k == 0:
                        point_limits[name] = [dmin, dmax]

                    # update limits
                    limits = point_limits[name]
                    if limits[0] > dmin:
                        limits[0] = dmin
                    if limits[1] < dmax:
                        limits[1] = dmax
                    point_limits[name] = limits

                # cell data limits
                for name, attribute_type in xdmf_info.cell_data.items():
                    # min and max from data
                    data = cell_data[name]
                    data = np.dstack(data)
                    data = data.squeeze()
                    dmin = data.min()
                    dmax = data.max()

                    limits: Tuple[float, float]
                    if k == 0:
                        cell_limits[name] = [dmin, dmax]

                    # update limits
                    limits = cell_limits[name]
                    if limits[0] > dmin:
                        limits[0] = dmin
                    if limits[1] < dmax:
                        limits[1] = dmax
                    cell_limits[name] = limits

        return DataLimits(
            limits={
                **point_limits,
                **cell_limits,
            }
        )


def xdmfs_from_febio(febio_dir: Path, xdmf_dir: Path, overwrite: bool = False) -> List[Path]:
    """Create the XDMF files from the raw FEBio simulation results.

    Processes all subdirectories which contain vtks. This allows
    to create the xdmfs for all simulations.
    """
    console.rule(title="XDMF from FEBio", align="left", style="white")
    # processes all folders with vtk

    # iterate over all directories and subdirectories to collect vtk directories
    vtk_dirs: List[Path] = []
    for d in sorted(list(febio_dir.rglob("*"))):
        vtk_paths = sorted(list(d.glob("*.vtk")))
        n_vtks = len(vtk_paths)
        if n_vtks > 0:
            console.print(f"{n_vtks} VTKs: {d}")
            vtk_dirs.append(d)

    if not vtk_dirs:
        logger.error(f"No directory with VTKs in '{febio_dir}'. "
                     f"Check if the 'febio_dir' is correct.")

    # process the VTKs in the directory
    xdmf_paths: List[Path] = []
    for vtk_dir in vtk_dirs:
        d_name = str(vtk_dir.relative_to(febio_dir))
        d_name = d_name.replace("/", "__")
        xdmf_path = xdmf_dir / f"{d_name}.xdmf"
        xdmf_paths.append(xdmf_path)
        vtks_to_xdmf(vtk_dir=vtk_dir, xdmf_path=xdmf_path, overwrite=overwrite)

    return xdmf_paths


def vtks_to_xdmf(vtk_dir: Path, xdmf_path: Path, overwrite: bool = False) -> None:
    """Convert VTK timesteps to XDMF time course.

    :param overwrite: forces overwrite of files, by default existing files are ignored.
    """
    console.rule(title=f"{xdmf_path}", style="white")

    vtk_paths = sorted(list(vtk_dir.glob("*.vtk")))
    mesh: meshio.Mesh = meshio.read(vtk_paths[0])

    if not xdmf_path.parent.exists():
        xdmf_path.parent.mkdir(parents=True)

    console.print(f"{vtk_dir} -> {xdmf_path}")
    if not overwrite and xdmf_path.exists():
        console.print(f"xdmf file exists: {xdmf_path}")
        return

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


def interpolate_xdmf(
    xdmf_in: Path, xdmf_out: Path, times_interpolate: np.ndarray
) -> None:
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
                raise ValueError(
                    f"Lower interpolation range outside of data: {times_interpolate[0]} < {times_data[0]}"
                )
            if times_interpolate[-1] > times_data[-1]:
                raise ValueError(
                    f"Upper interpolation range outside of data: {times_interpolate[-1]} > {times_data[-1]}"
                )

            lower_indices = np.zeros_like(times_interpolate, dtype=int)
            upper_indices = np.zeros_like(times_interpolate, dtype=int)
            for ki, ti in enumerate(times_interpolate):
                lower_indices[ki] = np.argwhere(times_data <= ti).flatten()[-1]
                upper_indices[ki] = np.argwhere(times_data >= ti).flatten()[0]

            # interpolate data for all data points
            for k in track(
                range(len(times_interpolate)), description="Interpolate data ..."
            ):
                t_interpolate = times_interpolate[k]
                idx_low = lower_indices[k]
                idx_up = upper_indices[k]

                # interpolate all numpy matrices
                t_low, point_data_low, cell_data_low = reader.read_data(idx_low)
                t_up, point_data_up, cell_data_up = reader.read_data(idx_up)
                f_up: float
                f_low: float
                if np.isclose(t_up, t_low):
                    f_low = 0.0
                    f_up = 1.0
                else:
                    f_low = (t_interpolate - t_low) / (t_up - t_low)
                    f_up = (t_up - t_interpolate) / (t_up - t_low)

                # interpolate point data
                point_data = {}
                for key in point_data_low:
                    point_data[key] = (f_low * point_data_low[key]) + (
                        f_up * point_data_up[key]
                    )

                # interpolate cell data
                cell_data = {}
                for key in cell_data_low:
                    # console.print(f"{cell_data_low[key]}")
                    # process the list
                    cell_data[key] = [
                        (f_low * cell_data_low[key][k]) + (f_up * cell_data_low[key][k])
                        for k in range(len(cell_data_low[key]))
                    ]

                # write interpolation point
                writer.write_data(
                    t_interpolate, point_data=point_data, cell_data=cell_data
                )

        # Fix incorrect *.h5 path
        # https://github.com/nschloe/meshio/pull/1358
        h5_path = xdmf_out.parent / f"{xdmf_out.stem}.h5"
        if h5_path.exists():
            os.remove(h5_path)
        shutil.move(f"{xdmf_out.stem}.h5", str(xdmf_out.parent))

        console.print(f"Interpolated data: {xdmf_out}")


if __name__ == "__main__":
    xdmf_path: Path = Path('/home/mkoenig/git/porous_media/data/spt_substrate_scan/sim_25.xdmf')
    xdmf_info: XDMFInformation = XDMFInformation.from_path(xdmf_path)
    console.print(xdmf_info)

    limits = data_limits_from_xdmf(xdmf_path)
    console.print(limits)
