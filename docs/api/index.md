# API reference

The reference is rendered from the installed source tree using mkdocstrings.
Use Python's `pathlib.Path` for filesystem arguments.

| Module area | Main entry points |
| --- | --- |
| [Data](data.md) | `vtks_to_xdmf`, `interpolate_xdmf`, `XDMFInfo`, `DataLimits` |
| [Mesh](mesh.md) | `ZonatedMesh`, `ZonationPatterns`, `mesh_to_vtk`, `mesh_to_xdmf` |
| [Visualization](visualization.md) | `DataLayer`, `VisualizationSettings`, image and video helpers |
| [FEBio](febio.md) | `create_feb_files` |

`porous_media.RESOURCES_DIR` points to resources distributed with the package.
Use your own output paths rather than relying on repository-relative constants
such as `DATA_DIR` and `RESULTS_DIR` in installed environments.

Research scripts under `porous_media.analyses` are intentionally outside this
reference. They can depend on unpublished input data and local paths.

The [user guides](../index.md#start-here) describe input conventions and
limitations that signatures alone cannot express.
