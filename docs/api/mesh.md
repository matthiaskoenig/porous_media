# Mesh

See [mesh zonation](../zonation.md) for boundary labels, volume requirements,
position direction, normalization, and the current single-cell-block assumption.

## Zonation

::: porous_media.mesh.mesh_zonation
    options:
      members:
        - ZonationPatterns
        - ZonatedMesh
        - visualize_zonation_patterns
        - example_mesh_zonation

## Serialization

These helpers serialize static meshes. Use the time-series conversion helpers
in the [data API](data.md) for simulation time courses.

::: porous_media.mesh.mesh_tools
    options:
      members:
        - mesh_to_xdmf
        - mesh_to_vtk
