"""Process timeseries data in XDMF.

https://www.xdmf.org/index.php/XDMF_Model_and_Format

"""
from __future__ import annotations

import json
import os
import shutil
from copy import deepcopy
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Any, Dict, List, Tuple

import meshio
import numpy as np
from dataclasses_json import dataclass_json
from meshio import CellBlock
from meshio._common import raw_from_cell_data
from meshio.xdmf.common import attribute_type
from rich.progress import track

from porous_media.console import console
from porous_media.log import get_logger


logger = get_logger(__name__)


class AttributeType(str, Enum):
    """Definition of variable type."""

    SCALAR = "Scalar"
    VECTOR = "Vector"
    TENSOR = "Tensor"
    MATRIX = "Matrix"
    OTHER = "Other"


@dataclass
class AttributeInfo:
    """Information on attribute on the mesh (cell or point data)."""

    name: str
    attribute_type: AttributeType
    shape: Tuple


@dataclass
class XDMFInfo:
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
    num_points: int
    num_cells: int
    points: List[Any]
    cells: List[Any]
    point_data: Dict[str, AttributeInfo]
    cell_data: Dict[str, AttributeInfo]

    @staticmethod
    def from_path(xdmf_path: Path) -> XDMFInfo:
        """Create XDMFInformation from xdmf path."""
        if not xdmf_path:
            raise IOError(f"xdmf_path does not exist: {xdmf_path}")

        with meshio.xdmf.TimeSeriesReader(xdmf_path) as reader:
            points, cells = reader.read_points_cells()

            # process most information based on initial data and points and cells
            tstart, point_data, cell_data = reader.read_data(0)

            point_data_info = {}
            for name, data in point_data.items():
                point_data_info[name] = AttributeInfo(
                    name=name,
                    attribute_type=attribute_type(data),
                    shape=np.dstack(data).shape,
                )

            cell_data_info = {}
            raw = raw_from_cell_data(cell_data)
            for name, data in raw.items():
                cell_data_info[name] = AttributeInfo(
                    name=name,
                    attribute_type=AttributeType(attribute_type(data)),
                    shape=np.dstack(data).shape,
                )

            # get end time
            tend, _, _ = reader.read_data(reader.num_steps - 1)

        cell_block: CellBlock = cells[0]
        xdmf_info = XDMFInfo(
            path=xdmf_path,
            tstart=tstart,
            tend=tend,
            num_steps=reader.num_steps,
            num_points=len(points),
            num_cells=len(cell_block),
            points=points,
            cells=cells,
            point_data=point_data_info,
            cell_data=cell_data_info,
        )
        return xdmf_info


@dataclass_json
@dataclass
class DataLimits:
    """Limits of numerical data."""

    limits: Dict[str, Tuple[float, float]]

    @classmethod
    def json_path_from_xdmf(cls, xdmf_path: Path) -> Path:
        """Calculate JSON path from xdmf path."""
        return xdmf_path.parent / f"{xdmf_path.stem}_limits.json"

    @classmethod
    def from_xdmf(cls, xdmf_path: Path, overwrite: bool = False) -> DataLimits:
        """Calculate data limits from XDMF timecourse..

        This takes some time because it requires iteration over the complete dataset.
        Should support subsets and operations such as merging of limits from multiple
        simulations.
        """
        # check existing Json
        json_path = cls.json_path_from_xdmf(xdmf_path)
        if not overwrite and json_path.exists():
            console.print(f"json file exists: {json_path}")
            with open(json_path, "r") as f_json:
                d = json.load(f_json)
                return DataLimits(**d)

        # read information
        xdmf_info = XDMFInfo.from_path(xdmf_path)

        point_limits: Dict[str, Tuple[float, float]] = {}
        cell_limits: Dict[str, Tuple[float, float]] = {}
        with meshio.xdmf.TimeSeriesReader(xdmf_path) as reader:
            _, _ = reader.read_points_cells()

            for k in track(
                range(reader.num_steps), description="Calculating data limits ..."
            ):
                t, point_data, cell_data = reader.read_data(k)

                # point data limits
                for name in xdmf_info.point_data:
                    # min and max from data
                    data = point_data[name]
                    data = np.dstack(data)
                    data = data.squeeze()
                    dmin = float(data.min())  # casting for JSON serialization
                    dmax = float(data.max())

                    limits: Tuple[float, float]
                    if k == 0:
                        point_limits[name] = (dmin, dmax)

                    # update limits
                    limits = point_limits[name]
                    point_limits[name] = (
                        dmin if limits[0] > dmin else limits[0],
                        dmax if limits[1] < dmax else limits[1],
                    )

                # cell data limits
                for name in xdmf_info.cell_data:
                    # min and max from data
                    data = cell_data[name]
                    data = np.dstack(data)
                    data = data.squeeze()
                    dmin = float(data.min())
                    dmax = float(data.max())

                    if k == 0:
                        cell_limits[name] = (dmin, dmax)

                    # update limits
                    limits = cell_limits[name]
                    cell_limits[name] = (
                        dmin if limits[0] > dmin else limits[0],
                        dmax if limits[1] < dmax else limits[1],
                    )

        data_limits = DataLimits(
            limits={
                **point_limits,
                **cell_limits,
            }
        )
        with open(json_path, "w") as f_json:
            djson: Dict = data_limits.to_dict()  # type: ignore
            json.dump(djson, fp=f_json, indent=2)

        console.print(f"json file created: {json_path}")
        return data_limits

    @staticmethod
    def merge_limits(data_limits: List[DataLimits]) -> DataLimits:
        """Merge the limits from multiple data limits."""

        # copy limits from simulation 0
        limits = deepcopy(data_limits[0].limits)
        for k, dlim in enumerate(data_limits):
            if k == 0:
                continue
            else:
                for name, lims_k in dlim.limits.items():
                    lims = limits[name]
                    limits[name] = (
                        lims_k[0] if lims_k[0] < lims[0] else lims[0],
                        lims_k[1] if lims_k[1] > lims[1] else lims[1],
                    )

        return DataLimits(limits=limits)


def xdmfs_from_directory(
    input_dir: Path, xdmf_dir: Path, overwrite: bool = False
) -> Dict[Path, Path]:
    """Create the XDMF files from a given input directory.

    E.g. this is the FEBIO output directory. This should be a directory with possible
    subdirectories which contain VTK files. All subdirectories are processed.

    Processes directory and all subdirectories which contain vtks. This allows to
    easily create all xdms for a given directory.
    XDMF file structure is analogue to the structure in the original directory

    :param input_dir: directory with VTKs, possible subdirectories
    :param xdmf_dir: directory in which the xdmfs are generated
    :returns: Dictionary {xdmf_path: subdirectory}
    """
    console.rule(title="XDMF from directory", align="left", style="white")
    # processes all folders with vtk

    # iterate over all directories and subdirectories to collect vtk directories
    vtk_dirs: List[Path] = []
    all_dirs = sorted([input_dir] + list(input_dir.rglob("*")))

    for d in all_dirs:
        vtk_paths = sorted(list(d.glob("*.vtk")))
        n_vtks = len(vtk_paths)
        if n_vtks > 0:
            console.print(f"{n_vtks} VTKs: {d}")
            vtk_dirs.append(d)

    if not vtk_dirs:
        logger.error(
            f"No directory with VTKs in '{input_dir}'. "
            f"Check if the 'febio_dir' is correct."
        )

    # process the VTKs in the directory
    xdmf_dict: Dict[Path, Path] = {}
    for vtk_dir in vtk_dirs:
        d_name = str(vtk_dir.relative_to(input_dir))
        d_name = d_name.replace("/", "__")
        xdmf_path = xdmf_dir / f"{d_name}.xdmf"
        xdmf_dict[xdmf_path] = vtk_dir
        vtks_to_xdmf(vtk_dir=vtk_dir, xdmf_path=xdmf_path, overwrite=overwrite)

    return xdmf_dict


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
        DataLimits.from_xdmf(xdmf_path=xdmf_path, overwrite=overwrite)
        return

    with meshio.xdmf.TimeSeriesWriter(xdmf_path, data_format="HDF") as writer:
        writer.write_points_cells(mesh.points, mesh.cells)

        for k in track(range(len(vtk_paths)), description="Processing VTKs ..."):
            vtk_path = vtk_paths[k]
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

    # Calculate limits
    DataLimits.from_xdmf(xdmf_path=xdmf_path, overwrite=overwrite)


def interpolate_xdmf(
    xdmf_in: Path,
    xdmf_out: Path,
    times_interpolate: np.ndarray,
    overwrite: bool = False,
) -> None:
    """Interpolate XDMF."""
    console.rule(title=f"Interpolate {xdmf_in}", style="white")

    if not overwrite and xdmf_out.exists():
        console.print(f"xdmf file exists: {xdmf_out}")
        DataLimits.from_xdmf(xdmf_path=xdmf_out, overwrite=overwrite)
        return

    with meshio.xdmf.TimeSeriesReader(xdmf_in) as reader:
        points, cells = reader.read_points_cells()

        with meshio.xdmf.TimeSeriesWriter(xdmf_out) as writer:
            writer.write_points_cells(points, cells)

            tnum = reader.num_steps
            times_data = np.zeros((tnum,))
            for k in track(range(tnum), description="Calculating timepoints ..."):
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
                range(len(times_interpolate)), description="Interpolating data ..."
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

        # Calculate limits
        DataLimits.from_xdmf(xdmf_path=xdmf_out, overwrite=overwrite)

        console.print(f"Interpolated data: {xdmf_out}")


if __name__ == "__main__":
    from porous_media import RESOURCES_DIR, RESULTS_DIR

    # Simple test that all variables are read
    # vtk_dir = RESOURCES_DIR / "vtk" / "vtk_single"
    vtk_dir = RESOURCES_DIR / "vtk" / "vtk_timecourse"
    xdmf_path = RESULTS_DIR / "vtk_test.xdmf"
    vtks_to_xdmf(vtk_dir, xdmf_path=xdmf_path, overwrite=True)
    xdmf_info: XDMFInfo = XDMFInfo.from_path(xdmf_path)
    console.print(xdmf_info)

    xdmf_paths: List[Path] = [
        Path(f"/home/mkoenig/git/porous_media/data/spt_substrate_scan/sim_{k}.xdmf")
        for k in range(21, 26)
    ]

    xdmf_path = Path(
        "/home/mkoenig/git/porous_media/data/spt_substrate_scan/sim_25.xdmf"
    )
    xdmf_info = XDMFInfo.from_path(xdmf_path)
    console.print(xdmf_info)

    # limits of individual simulations
    all_limits: List[DataLimits] = []
    for xdmf_path in xdmf_paths:
        limits = DataLimits.from_xdmf(xdmf_path=xdmf_path, overwrite=False)
        all_limits.append(limits)

    # merge data limits
    data_limits = DataLimits.merge_limits(all_limits)
    console.print(data_limits)
