"""Process timeseries data in XDMF.

https://www.xdmf.org/index.php/XDMF_Model_and_Format

"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path
from typing import cast
from xml.etree import ElementTree as ET

import h5py
import meshio
import numpy as np
from dataclasses_json import dataclass_json
from meshio._common import raw_from_cell_data
from meshio._exceptions import ReadError, WriteError
from meshio.xdmf.common import attribute_type, numpy_to_xdmf_dtype
from rich.progress import track

from porous_media.console import console
from porous_media.data.netcdf import HDF_COMPRESSION
from porous_media.log import get_logger

logger = get_logger(__name__)

CHUNK_TARGET_BYTES = 1 << 20
"""Target size of one HDF chunk of a stacked attribute. A chunk spanning
several time steps is what lets the compressor exploit how little changes
between them; 1 MiB measured 3.5x smaller than one dataset per step while
keeping a single-step read near 4 ms."""

MAX_CHUNK_STEPS = 512
"""Upper bound on the time steps in one chunk, for attributes small enough that
:data:`CHUNK_TARGET_BYTES` would span more steps than a series contains."""


class XDMFTimeSeriesWriter(meshio.xdmf.TimeSeriesWriter):
    """Write time series with compressed HDF data beside XDMF, in place.

    This preserves meshio's ``TimeSeriesWriter`` context manager API and writes
    relative HDF references so that the XDMF and HDF files can be moved together.

    Each HDF dataset is deflated with the byte-shuffle filter. This is a
    standard HDF5 filter; readers need nothing special and the stored values and
    dtypes are unchanged.

    With ``share_identical_arrays``, an array identical to one already written
    shares its storage through an HDF link instead of being stored twice. A
    fixed mesh repeats its time-invariant fields at every step, so sharing keeps
    a long history proportional to the data that actually changes. It is not the
    default because it makes the sidecar unsafe to edit in place: writing to one
    dataset then changes every array that shares its storage.

    With ``stack_time_steps``, each attribute becomes one growing dataset
    spanning the whole series, and every step references its slice through an
    XDMF HyperSlab. A dataset written one step at a time gives the compressor
    only a single frame to work with, while a chunk spanning consecutive steps
    exposes how little changes between them: on reconstructed scalar fields this
    measured 3.5x smaller than one dataset per step. The chunk spans as many
    steps as fit in :data:`CHUNK_TARGET_BYTES`, which keeps a single-step read
    to a few milliseconds.

    meshio cannot read HyperSlab references, so a stacked series must be read
    with :class:`XDMFTimeSeriesReader`. ParaView reads them directly.
    """

    def __init__(
        self,
        filename: str | Path,
        data_format: str = "HDF",
        *,
        share_identical_arrays: bool = False,
        stack_time_steps: bool = False,
    ) -> None:
        """Record how arrays of the series are laid out in the HDF sidecar."""
        super().__init__(filename, data_format=data_format)
        self.share_identical_arrays = share_identical_arrays
        self.stack_time_steps = stack_time_steps

    def __enter__(self) -> XDMFTimeSeriesWriter:
        """Open the HDF sidecar when using HDF output."""
        # meshio writes its HDF file into the working directory by default.
        # Keep the relative filename in the XML, but open it at its actual location.
        self._datasets: dict[tuple[str, tuple[int, ...], bytes], str] = {}
        self._stacks: dict[tuple[str, str], h5py.Dataset] = {}
        self._stack_paths: set[str] = set()
        self._stack_sources: list[tuple[ET.Element, h5py.Dataset]] = []
        if self.data_format == "HDF":
            self.h5_filename = str(self.filename.with_suffix(".h5"))
            self.h5_file = h5py.File(self.h5_filename, "w")
        return self

    def numpy_to_xml_string(self, data: np.ndarray) -> str:
        """Write one array compressed, optionally reusing an identical array."""
        if self.data_format != "HDF":
            return super().numpy_to_xml_string(data)
        data = np.ascontiguousarray(data)
        name = f"data{self.data_counter}"
        self.data_counter += 1
        previous = None
        if self.share_identical_arrays:
            key = (
                data.dtype.str,
                data.shape,
                hashlib.blake2b(data, digest_size=32).digest(),
            )
            previous = self._datasets.get(key)
            if previous is None:
                self._datasets[key] = name
        if previous is None:
            self.h5_file.create_dataset(name, data=data, **HDF_COMPRESSION)
        else:
            self.h5_file[name] = self.h5_file[previous]
        return f"{Path(self.h5_filename).name}:/{name}"

    def _stacked_dataset(
        self, center: str, name: str, data: np.ndarray
    ) -> h5py.Dataset:
        """Return the growing dataset of one attribute, extended by this step."""
        key = (center, name)
        dataset = self._stacks.get(key)
        if dataset is None:
            # Keep a readable, collision-free path for the attribute.
            path = f"{center}/{name.replace('/', '_')}"
            while path in self._stack_paths:
                path = f"{path}_"
            self._stack_paths.add(path)
            # The series length is unknown while writing, so the span is capped:
            # a chunk far longer than the series would be mostly empty padding.
            span = max(
                1, min(MAX_CHUNK_STEPS, CHUNK_TARGET_BYTES // max(data.nbytes, 1))
            )
            dataset = self.h5_file.create_dataset(
                path,
                shape=(0, *data.shape),
                maxshape=(None, *data.shape),
                dtype=data.dtype,
                chunks=(span, *data.shape),
                **HDF_COMPRESSION,
            )
            self._stacks[key] = dataset
        if dataset.shape[1:] != data.shape or dataset.dtype != data.dtype:
            raise WriteError(
                f"Attribute '{name}' changes shape or type between time steps"
            )
        dataset.resize(dataset.shape[0] + 1, axis=0)
        dataset[-1] = data
        return dataset

    def _stacked_attribute(
        self, grid: ET.Element, name: str, data: np.ndarray, center: str
    ) -> None:
        """Reference this step's slice of the attribute's dataset as a HyperSlab."""
        attribute = ET.SubElement(
            grid,
            "Attribute",
            Name=name,
            AttributeType=attribute_type(data),
            Center=center,
        )
        dataset = self._stacked_dataset(center, name, data)
        step = dataset.shape[0] - 1
        rank = len(dataset.shape)
        item = ET.SubElement(
            attribute,
            "DataItem",
            ItemType="HyperSlab",
            Dimensions=" ".join(str(size) for size in data.shape),
            Type="HyperSlab",
        )
        selection = ET.SubElement(
            item, "DataItem", Dimensions=f"3 {rank}", Format="XML", DataType="Int"
        )
        rows = (
            [step, *([0] * (rank - 1))],  # start
            [1] * rank,  # stride
            [1, *data.shape],  # count
        )
        selection.text = "\n".join(" ".join(str(v) for v in row) for row in rows)
        data_type, precision = numpy_to_xdmf_dtype[data.dtype.name]
        source = ET.SubElement(
            item,
            "DataItem",
            DataType=data_type,
            # The series is still growing; the final extent is set on close.
            Dimensions=" ".join(str(size) for size in dataset.shape),
            Format="HDF",
            Precision=precision,
        )
        source.text = f"{Path(self.h5_filename).name}:/{dataset.name.lstrip('/')}"
        self._stack_sources.append((source, dataset))

    def point_data(self, point_data: dict[str, np.ndarray], grid: ET.Element) -> None:
        """Write node attributes, stacking them when the series is stacked."""
        if not self._stacking():
            super().point_data(point_data, grid)
            return
        for name, data in point_data.items():
            self._stacked_attribute(grid, name, np.asarray(data), "Node")

    def cell_data(
        self, cell_data: dict[str, list[np.ndarray]], grid: ET.Element
    ) -> None:
        """Write cell attributes, stacking them when the series is stacked."""
        if not self._stacking():
            super().cell_data(cell_data, grid)
            return
        for name, data in raw_from_cell_data(cell_data).items():
            self._stacked_attribute(grid, name, np.asarray(data), "Cell")

    def _stacking(self) -> bool:
        """Report whether attributes are written as slices of one dataset."""
        return self.stack_time_steps and self.data_format == "HDF"

    def __exit__(self, *args: object) -> None:
        """Close the HDF sidecar even if serializing XDMF metadata fails."""
        try:
            # Every step is written, so each stacked reference can state the
            # final extent of its dataset.
            for source, dataset in self._stack_sources:
                source.set("Dimensions", " ".join(str(s) for s in dataset.shape))
            super().__exit__(*args)
        finally:
            if self.data_format == "HDF":
                self.h5_file.close()


class XDMFTimeSeriesReader(meshio.xdmf.TimeSeriesReader):
    """Read a time series, resolving the HyperSlab references meshio cannot.

    :class:`XDMFTimeSeriesWriter` can store each attribute as one dataset
    spanning the whole series, with every step referencing its slice through an
    XDMF HyperSlab. meshio has no HyperSlab support, so reading such a series
    needs this subclass. A series without HyperSlab references is read exactly
    as meshio reads it, which makes this a drop-in replacement.

    Only the requested slice is read, so a step costs one chunk rather than the
    whole history.
    """

    def _read_data_item(self, data_item: ET.Element) -> np.ndarray:
        """Resolve one data item, slicing a stacked dataset when referenced."""
        if data_item.get("ItemType") != "HyperSlab":
            return super()._read_data_item(data_item)
        children = list(data_item)
        if len(children) != 2:
            raise ReadError("A HyperSlab needs a selection and a source data item")
        selection_item, source_item = children
        values = [int(value) for value in (selection_item.text or "").split()]
        if len(values) % 3:
            raise ReadError("A HyperSlab selection needs start, stride and count")
        rank = len(values) // 3
        start, stride, count = (
            values[:rank],
            values[rank : 2 * rank],
            values[2 * rank :],
        )
        source = self._hdf_dataset(source_item)
        if source is None:
            # An inline or binary source is small enough to resolve whole.
            whole = super()._read_data_item(source_item)
            selected = whole[
                tuple(
                    slice(begin, begin + size * step, step)
                    for begin, step, size in zip(start, stride, count, strict=True)
                )
            ]
        else:
            if len(source.shape) != rank:
                raise ReadError("A HyperSlab selection must match the source rank")
            selected = source[
                tuple(
                    slice(begin, begin + size * step, step)
                    for begin, step, size in zip(start, stride, count, strict=True)
                )
            ]
        dims = [int(d) for d in (data_item.get("Dimensions") or "").split()]
        return np.asarray(selected).reshape(dims) if dims else np.asarray(selected)

    def _hdf_dataset(self, data_item: ET.Element) -> h5py.Dataset | None:
        """Return the HDF dataset a data item names, or None if it is not HDF."""
        if data_item.attrib.get("Format") != "HDF" or not data_item.text:
            return None
        filename, _, h5path = data_item.text.strip().partition(":")
        if not h5path.startswith("/"):
            raise ReadError(f"HDF reference must be absolute: {data_item.text!r}")
        # The HDF path is relative to the XDMF file, so the pair stays movable.
        full = self.filename.resolve().parent / filename
        handle = self.hdf5_files.get(full)
        if handle is None:
            handle = h5py.File(full, "r")
            self.hdf5_files[full] = handle
        node: object = handle
        for key in h5path[1:].split("/"):
            node = cast("h5py.Group", node)[key]
        if not isinstance(node, h5py.Dataset):
            raise ReadError(f"HDF reference is not a dataset: {data_item.text!r}")
        return node


class AttributeType(StrEnum):
    """Definition of variable type."""

    SCALAR = "Scalar"
    VECTOR = "Vector"
    TENSOR = "Tensor"
    TENSOR6 = "Tensor6"
    MATRIX = "Matrix"
    OTHER = "Other"


@dataclass
class AttributeInfo:
    """Information on attribute on the mesh (cell or point data)."""

    name: str
    attribute_type: AttributeType
    shape: tuple[int, ...]


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
    points: np.ndarray
    cells: list[meshio.CellBlock]
    point_data: dict[str, AttributeInfo]
    cell_data: dict[str, AttributeInfo]

    @staticmethod
    def from_path(xdmf_path: Path) -> XDMFInfo:
        """Create XDMFInformation from xdmf path."""
        if not xdmf_path.is_file():
            raise FileNotFoundError(f"XDMF path does not exist: {xdmf_path}")

        with XDMFTimeSeriesReader(xdmf_path) as reader:
            points, cells = reader.read_points_cells()
            if reader.num_steps == 0:
                raise ValueError("XDMF time series must contain at least one time step")

            # process most information based on initial data and points and cells
            tstart, point_data, cell_data = reader.read_data(0)

            point_data_info = {}
            for name, data in point_data.items():
                point_data_info[name] = AttributeInfo(
                    name=name,
                    attribute_type=AttributeType(attribute_type(data)),
                    shape=data.shape,
                )

            cell_data_info = {}
            raw = raw_from_cell_data(cell_data)
            for name, data in raw.items():
                cell_data_info[name] = AttributeInfo(
                    name=name,
                    attribute_type=AttributeType(attribute_type(data)),
                    shape=data.shape,
                )

            # get end time
            tend, _, _ = reader.read_data(reader.num_steps - 1)

        return XDMFInfo(
            path=xdmf_path,
            tstart=tstart,
            tend=tend,
            num_steps=reader.num_steps,
            num_points=len(points),
            num_cells=sum(len(cell_block) for cell_block in cells),
            points=points,
            cells=cells,
            point_data=point_data_info,
            cell_data=cell_data_info,
        )


@dataclass_json
@dataclass
class DataLimits:
    """Limits of numerical data."""

    limits: dict[str, tuple[float, float]]

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
            with json_path.open(encoding="utf-8") as f_json:
                d = json.load(f_json)
                return DataLimits(
                    limits={name: tuple(values) for name, values in d["limits"].items()}
                )

        # read information
        xdmf_info = XDMFInfo.from_path(xdmf_path)

        point_limits: dict[str, tuple[float, float]] = {}
        cell_limits: dict[str, tuple[float, float]] = {}
        with XDMFTimeSeriesReader(xdmf_path) as reader:
            _, _ = reader.read_points_cells()

            for k in track(
                range(reader.num_steps), description="Calculating data limits ..."
            ):
                _t, point_data, cell_data = reader.read_data(k)

                # point data limits
                for name in xdmf_info.point_data:
                    # min and max from data
                    data = point_data[name]
                    data = np.asarray(data)
                    data = data.squeeze()
                    dmin = float(data.min())  # casting for JSON serialization
                    dmax = float(data.max())

                    limits: tuple[float, float]
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
                    data = np.concatenate(cell_data[name])
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
        with json_path.open("w", encoding="utf-8") as f_json:
            json.dump({"limits": data_limits.limits}, fp=f_json, indent=2)

        console.print(f"json file created: {json_path}")
        return data_limits

    @staticmethod
    def merge_limits(data_limits: list[DataLimits]) -> DataLimits:
        """Merge the limits from multiple data limits."""
        if not data_limits:
            raise ValueError("data_limits must contain at least one set of limits")
        limits: dict[str, tuple[float, float]] = {}
        for dlim in data_limits:
            for name, lims_k in dlim.limits.items():
                lims = limits.get(name, lims_k)
                limits[name] = (
                    lims_k[0] if lims_k[0] < lims[0] else lims[0],
                    lims_k[1] if lims_k[1] > lims[1] else lims[1],
                )

        return DataLimits(limits=limits)


def xdmfs_from_directory(
    input_dir: Path, xdmf_dir: Path, overwrite: bool = False
) -> dict[Path, Path]:
    """Create the XDMF files from a given input directory.

    E.g. this is the FEBIO output directory. This should be a directory with possible
    subdirectories which contain VTK files. All subdirectories are processed.

    Processes directory and all subdirectories which contain vtks. This allows to
    easily create all xdms for a given directory.
    XDMF file structure is analogue to the structure in the original directory

    Args:
        input_dir: Directory containing VTK files, including subdirectories.
        xdmf_dir: Directory in which the XDMF files are generated.
        overwrite: Whether to replace existing output files.

    Returns:
        A mapping from each XDMF output path to its VTK source directory.
    """
    if not input_dir.is_dir():
        raise NotADirectoryError(f"VTK input directory does not exist: {input_dir}")
    vtk_dirs = sorted(
        {path.parent for path in input_dir.rglob("*.vtk") if path.is_file()}
    )
    if not vtk_dirs:
        raise ValueError(f"No VTK files found in {input_dir!s}")
    outputs: dict[Path, Path] = {}
    for vtk_dir in vtk_dirs:
        relative = vtk_dir.relative_to(input_dir)
        name = "__".join(relative.parts) if relative.parts else input_dir.resolve().name
        output = xdmf_dir / f"{name}.xdmf"
        if output in outputs:
            raise ValueError(f"VTK directories map to the same output: {output}")
        outputs[output] = vtk_dir
    for output, vtk_dir in outputs.items():
        vtks_to_xdmf(vtk_dir, output, overwrite=overwrite)
    return outputs


def vtks_to_xdmf(vtk_dir: Path, xdmf_path: Path, overwrite: bool = False) -> None:
    """Convert VTK timesteps to XDMF time course.

    Args:
        vtk_dir: Directory containing VTK files ordered by timestep.
        xdmf_path: Path to the XDMF output file.
        overwrite: Whether to replace existing output files.
    """
    console.rule(title=f"{xdmf_path}", style="white")

    console.print(f"{vtk_dir} -> {xdmf_path}")
    if not overwrite and xdmf_path.exists():
        console.print(f"xdmf file exists: {xdmf_path}")
        DataLimits.from_xdmf(xdmf_path=xdmf_path, overwrite=overwrite)
        return

    vtk_paths = sorted(vtk_dir.glob("*.vtk"))
    if not vtk_paths:
        raise ValueError(f"No VTK files found in '{vtk_dir}'")
    if xdmf_path.is_symlink() or xdmf_path.with_suffix(".h5").is_symlink():
        raise ValueError("XDMF outputs must not follow symbolic links")
    if not overwrite and xdmf_path.with_suffix(".h5").exists():
        raise FileExistsError(
            f"HDF sidecar already exists: {xdmf_path.with_suffix('.h5')}"
        )
    mesh: meshio.Mesh = meshio.read(vtk_paths[0])
    xdmf_path.parent.mkdir(parents=True, exist_ok=True)

    # A converted series keeps the plain per-step layout every XDMF reader
    # understands; the fields that never change share one array instead.
    with XDMFTimeSeriesWriter(
        xdmf_path, data_format="HDF", share_identical_arrays=True
    ) as writer:
        writer.write_points_cells(
            mesh.points, [(cell.type, cell.data) for cell in mesh.cells]
        )

        for k in track(range(len(vtk_paths)), description="Processing VTKs ..."):
            vtk_path = vtk_paths[k]
            # read timepoint
            with vtk_path.open(encoding="utf-8") as f_vtk:
                f_vtk.readline()  # skip first line
                line = f_vtk.readline()
                t = float(line.split()[-1])

            # read mesh
            mesh = meshio.read(vtk_path)
            writer.write_data(t, point_data=mesh.point_data, cell_data=mesh.cell_data)

    # Calculate limits
    DataLimits.from_xdmf(xdmf_path=xdmf_path, overwrite=True)


def interpolate_xdmf(
    xdmf_in: Path,
    xdmf_out: Path,
    times_interpolate: np.ndarray,
    overwrite: bool = False,
) -> None:
    """Linearly interpolate point and cell data at finite times within the input range.

    The input time series must have strictly increasing, finite timepoints.
    Requested times must form a nonempty one-dimensional array. Their order and
    repeated values are preserved. Output directories are created as needed.
    """
    console.rule(title=f"Interpolate {xdmf_in}", style="white")

    if not overwrite and xdmf_out.exists():
        console.print(f"xdmf file exists: {xdmf_out}")
        DataLimits.from_xdmf(xdmf_path=xdmf_out, overwrite=overwrite)
        return

    times_interpolate = np.asarray(times_interpolate, dtype=float)
    if times_interpolate.ndim != 1 or times_interpolate.size == 0:
        raise ValueError("Interpolation times must be a nonempty one-dimensional array")
    if not np.isfinite(times_interpolate).all():
        raise ValueError("Interpolation times must be finite")
    if xdmf_in.resolve() == xdmf_out.resolve():
        raise ValueError("Input and output XDMF paths must differ")
    if xdmf_in.with_suffix(".h5").resolve() == xdmf_out.with_suffix(".h5").resolve():
        raise ValueError("Input and output HDF sidecar paths must differ")
    if xdmf_out.is_symlink() or xdmf_out.with_suffix(".h5").is_symlink():
        raise ValueError("XDMF outputs must not follow symbolic links")
    if not overwrite and xdmf_out.with_suffix(".h5").exists():
        raise FileExistsError(
            f"HDF sidecar already exists: {xdmf_out.with_suffix('.h5')}"
        )

    with XDMFTimeSeriesReader(xdmf_in) as reader:
        points, cells = reader.read_points_cells()

        times_data = np.zeros(reader.num_steps)
        for k in track(
            range(reader.num_steps), description="Calculating timepoints ..."
        ):
            t, _, _ = reader.read_data(k)
            times_data[k] = t

        if (
            times_data.size == 0
            or not np.isfinite(times_data).all()
            or not (np.diff(times_data) > 0).all()
        ):
            raise ValueError("Input timepoints must be finite and strictly increasing")
        if times_interpolate.min() < times_data[0]:
            raise ValueError(
                f"Lower interpolation range outside of data: {times_interpolate.min()} < {times_data[0]}"
            )
        if times_interpolate.max() > times_data[-1]:
            raise ValueError(
                f"Upper interpolation range outside of data: {times_interpolate.max()} > {times_data[-1]}"
            )

        xdmf_out.parent.mkdir(parents=True, exist_ok=True)
        with XDMFTimeSeriesWriter(xdmf_out, share_identical_arrays=True) as writer:
            writer.write_points_cells(points, cells)

            # interpolate data for all data points
            for k in track(
                range(len(times_interpolate)), description="Interpolating data ..."
            ):
                t_interpolate = times_interpolate[k]
                idx_up = int(np.searchsorted(times_data, t_interpolate))
                idx_low = idx_up if times_data[idx_up] == t_interpolate else idx_up - 1

                # interpolate all numpy matrices
                t_low, point_data_low, cell_data_low = reader.read_data(idx_low)
                t_up, point_data_up, cell_data_up = reader.read_data(idx_up)
                f_up: float
                f_low: float
                if idx_up == idx_low:
                    f_low = 1.0
                    f_up = 0.0
                else:
                    f_up = (t_interpolate - t_low) / (t_up - t_low)
                    f_low = 1.0 - f_up

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
                        (f_low * cell_data_low[key][k]) + (f_up * cell_data_up[key][k])
                        for k in range(len(cell_data_low[key]))
                    ]

                # write interpolation point
                writer.write_data(
                    t_interpolate, point_data=point_data, cell_data=cell_data
                )

        # Calculate limits
        DataLimits.from_xdmf(xdmf_path=xdmf_out, overwrite=True)

        console.print(f"Interpolated data: {xdmf_out}")
