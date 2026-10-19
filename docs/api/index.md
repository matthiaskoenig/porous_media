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
Pass explicit data and output paths. The library no longer exports
repository-relative `BASE_DIR`, `DATA_DIR`, `RESULTS_DIR`, or `EXAMPLE_*` constants.
Research path configuration lives in the checkout's `research` namespace.

Research scripts under `research/` are intentionally outside this
reference and release artifacts. They can depend on unpublished inputs and local
paths. See [research setup and migration](../research.md).

The [user guides](../index.md#start-here) describe input conventions and
limitations that signatures alone cannot express.
