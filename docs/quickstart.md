# Quickstart

This example uses files distributed with the package. It needs no external
simulation data, FEBio installation, graphics backend, or FFmpeg.

## Convert bundled VTK results

Save the following as `quickstart.py`:

```python
from pathlib import Path

from porous_media import RESOURCES_DIR
from porous_media.data.xdmf_tools import XDMFInfo, vtks_to_xdmf

output_dir = Path("results") / "quickstart"
output = output_dir / "example.xdmf"
vtks_to_xdmf(
    vtk_dir=RESOURCES_DIR / "vtk" / "vtk_timecourse",
    xdmf_path=output,
)

info = XDMFInfo.from_path(output)
print(f"{info.num_steps} steps: {info.tstart:g} to {info.tend:g}")
print("Point fields:", sorted(info.point_data))
print("Cell fields:", sorted(info.cell_data))
assert info.num_steps == 3
assert info.tstart == 0.0
assert info.tend == 220.0
```

Run it in the environment where you installed `porous_media`:

```bash
python quickstart.py
```

In a development checkout, use `uv run --no-sync python quickstart.py` after
the [environment setup](development.md#environment).

The summary includes `3 steps: 0 to 220`.
The output directory contains:

```text
results/quickstart/
├── example.xdmf
├── example.h5
└── example_limits.json
```

Keep the XDMF and HDF5 files together. XDMF stores metadata and references the
numerical arrays in its HDF5 companion. The JSON file caches field limits.
By default, conversion reuses an existing output. Pass `overwrite=True` when
you intentionally regenerate it.

## Read one time step

Append this block to the script:

```python
import meshio

with meshio.xdmf.TimeSeriesReader(output) as reader:
    points, cells = reader.read_points_cells()
    time, point_data, cell_data = reader.read_data(0)

print(f"First time: {time:g}")
print("Pressure shape:", point_data["effective_fluid_pressure_TPM"].shape)
print("Cell block count:", len(cells))
assert len(point_data["effective_fluid_pressure_TPM"]) == len(points)
```

`point_data` maps each field to a NumPy array indexed by point. `cell_data`
maps each field to a list of arrays, one per mesh cell block.

Continue with [VTK and XDMF data](xdmf.md) for interpolation and caching,
[mesh zonation](zonation.md) for spatial patterns, or
[visualization](visualization.md) for image export.
