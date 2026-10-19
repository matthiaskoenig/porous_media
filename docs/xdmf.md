# VTK and XDMF data

The data helpers work with meshio XDMF time series that share a fixed geometry
and topology. Point fields and cell fields remain separate.

## Convert a VTK time series

```python
from pathlib import Path

from porous_media import RESOURCES_DIR
from porous_media.data.xdmf_tools import vtks_to_xdmf

output = Path("results") / "simulation.xdmf"
vtks_to_xdmf(
    vtk_dir=RESOURCES_DIR / "vtk" / "vtk_timecourse",
    xdmf_path=output,
    overwrite=True,
)
```

The converter reads legacy `*.vtk` files in sorted filename order. It reads the
time from the last whitespace-separated token of each file's second line, as
used by the bundled FEBio exports. Use filenames whose lexical order matches
time order. Each file must describe the same mesh and compatible fields.
This function is not a converter for arbitrary VTK time metadata or moving meshes.

The output parent directory is created when needed. Conversion writes:

- `simulation.xdmf`: the time series metadata;
- `simulation.h5`: numerical arrays, beside the XDMF file;
- `simulation_limits.json`: cached limits across all time steps.

For multiple simulation directories, `xdmfs_from_directory(input_dir, xdmf_dir)`
finds directories containing `*.vtk` files. Output names flatten nested
directory names using `__`; its return value maps output paths to input directories.
Files directly inside the input directory use that directory's name. Name
collisions are rejected before conversion, including on Windows. An empty or
missing source directory raises an error.

With `overwrite=False`, conversion and interpolation preserve an existing XDMF
output and refuse to replace an orphaned HDF sidecar. New outputs cannot be
symbolic links. Regenerating a time series also regenerates its limits cache.

## Inspect metadata and limits

Continue using `output` from the conversion example:

```python
from porous_media.data.xdmf_tools import DataLimits, XDMFInfo

info = XDMFInfo.from_path(output)
print(info.num_steps, info.tstart, info.tend)
print(sorted(info.point_data))
print(sorted(info.cell_data))

limits = DataLimits.from_xdmf(output)
print(limits.limits["effective_fluid_pressure_TPM"])
```

`XDMFInfo` inspects the mesh and field metadata from the first time step and
reads the final time. `DataLimits` scans all time steps on its first call,
then reuses the adjacent JSON cache.
Attribute shapes describe the original arrays, without added stacking dimensions.
Use `DataLimits.from_xdmf(output, overwrite=True)` after modifying source data.
Field names must be distinct between point and cell data when sharing one
limits dictionary. For vector and tensor fields, limits are component extrema,
not vector magnitudes.

## Interpolate in time

```python
import numpy as np

from porous_media.data.xdmf_tools import interpolate_xdmf

interpolate_xdmf(
    xdmf_in=output,
    xdmf_out=output.with_name("simulation_interpolated.xdmf"),
    times_interpolate=np.array([0.0, 40.0, 80.0, 150.0, 220.0]),
    overwrite=True,
)
```

Interpolation is linear for numerical point and cell arrays and preserves
the mesh. Supply a nonempty array of requested times within the source time
range. Requested order and repeated values are preserved; source time steps
must be strictly increasing.
The function does not extrapolate, remesh, or infer physically meaningful
interpolation for categorical fields. Source endpoints are valid requested times.

## Scalar datasets

`porous_media.data.xdmf_calculations.mesh_datasets_from_xdmf` returns a
`(cell_dataset, point_dataset)` pair of xarray datasets and caches them as
`<stem>_cells.nc` and `<stem>_points.nc`. This helper selects scalar fields;
vectors and tensors are skipped. Cached datasets do not retain mesh geometry.
Scalars from all cell blocks are concatenated in mesh order, including blocks
with different sizes. Dimensions are `(time, cell)` and `(time, point)`, with
explicit coordinates. A single cell retains its dimension. Scalar names and
sizes must remain consistent across time steps; changing numeric types are
promoted without truncating fractional values.

Returned datasets are fully loaded and do not keep cache files open. Set
`overwrite=True` after changing the source. Directory conversion accepts a
`max_workers` limit, bounds the default by available CPUs and file count, and
does not return large datasets through its process pool. Use `max_workers=1`
for serial execution. Parallel workers use spawn to avoid forking threaded native
libraries. Put script calls that start workers under `if __name__ == "__main__":`;
use serial execution in interactive sessions. For lazy access to large existing
caches, use `xarray.open_dataset` inside an explicit context manager instead of
this fully loaded helper. Use meshio directly for geometry, vectors, tensors, or
individual cell blocks.

See the [data API reference](api/data.md) for signatures.
