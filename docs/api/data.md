# Data

See [VTK and XDMF data](../xdmf.md) for fixed-mesh assumptions, VTK time
metadata, interpolation, and cache behavior.

## XDMF tools

::: porous_media.data.xdmf_tools
    options:
      members:
        - AttributeType
        - AttributeInfo
        - XDMFInfo
        - DataLimits
        - xdmfs_from_directory
        - vtks_to_xdmf
        - interpolate_xdmf

## Scalar datasets

These helpers select scalar fields and do not preserve mesh geometry.

::: porous_media.data.xdmf_calculations
    options:
      members:
        - mesh_datasets_from_xdmf
        - mesh_datasets_from_xdmf_dir
