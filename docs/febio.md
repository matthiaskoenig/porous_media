# FEBio workflows

`porous_media` prepares FEBio input templates and processes exported results.
Install the solver, FEBio Studio, and any required research plugins separately.
The package does not execute simulations or decode `.xplt` databases directly.

## Render input templates

Use `create_feb_files` to render one file per simulation context:

```python
from pathlib import Path

from porous_media.febio.febio_manipulation import create_feb_files

work_dir = Path("results") / "febio"
work_dir.mkdir(parents=True, exist_ok=True)
template = work_dir / "model.feb.template"
template.write_text(
    '<febio_spec version="4.0">\n'
    "  <!-- Example parameter: {{ boundary_pressure }} -->\n"
    "</febio_spec>\n",
    encoding="utf-8",
)
create_feb_files(
    feb_template=template,
    feb_data={
        "low": {"boundary_pressure": 133},
        "high": {"boundary_pressure": 266},
    },
    output_dir=work_dir / "inputs",
)
```

This demonstrates template rendering and creates `model_low.feb` and
`model_high.feb`. The minimal template is not a complete simulation model.
Replace it with a validated model containing the geometry, materials, boundary
conditions, units, and plugins required for your study.

Template variables use Jinja2 syntax, for example `{{ boundary_pressure }}`.
The mapping keys identify simulations and become part of output filenames.
The filename prefix is the template name before its first dot.
Simulation identifiers must be single portable filename components, without
directory separators, traversal, or Windows device names. This prevents a context
key from writing outside the selected output directory. Generated output files
cannot be symbolic links. Existing files are preserved by default; pass
`overwrite=True` to replace regular output files explicitly.

Missing template variables now raise an error. Context values are XML-escaped,
and an immutable Jinja2 sandbox blocks private Python attributes. Templates and
context objects must still be trusted and reviewed: sandboxing does not impose
CPU or memory limits or prevent side effects in supplied Python objects. Avoid
passing callables or objects with privileged methods. The helper does not validate
the resulting FEBio model; validate it with the solver before use.

When migrating existing templates, supply every required variable and remove
assumptions that unescaped strings can insert XML. Use Jinja2's `safe` filter only
for reviewed XML fragments. Parameter contexts are no longer printed to output.

## Run the solver and export results

1. Open and validate the rendered model using your FEBio installation.
2. Run the model using the solver version and plugins required by that model.
3. Export the desired time steps and fields from FEBio Studio as legacy VTK
   files, keeping one simulation in each directory.
4. Convert the exported files using `vtks_to_xdmf` and inspect their time range
   and fields before analysis.

Use the [official FEBio resources](https://febio.org/) for solver installation,
current commands, model specifications, and export support. Solver binaries and
plugin compatibility depend on the installation.

The converter expects the final token of each legacy VTK file's second line to
contain the numerical time. Keep filenames ordered by time and preserve a fixed
mesh. See [VTK and XDMF data](xdmf.md) for the precise input assumptions and
[visualization](visualization.md) for image and video export.

## Historical research environments

The [historical run notes](run_febio_simulations.md) and
[RoadRunner/FEBio build notes](build_roadrunner_febio.md) preserve an older
workstation and cluster environment. They include pinned upstream revisions,
local paths, source changes, and private research plugins. They are not required
to install or use the Python package and are not current general installation
instructions.

See the [FEBio API reference](api/febio.md) and
[research archive](research.md) for related material.
