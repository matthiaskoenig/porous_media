"""Helper functions for calculating information on the Timecourse XDMF.

Includes things such as calculation of sum, ...
How to best handle this?
Provide some simple helpers for common operations, e.g.
- integral over the mesh (e.g. total necrosis, fat)
- relative fraction of the mesh (with max values or one)
- densities (normalization on volumes)

This should work on cell_data and point_data of the mesh

"""
from porous_media.console import console

if __name__ == "__main__":
    console.rule(title="XDMF calculations", style="white")

    # FIXME: use the example of one of the zonation patterns with 10 timepoints
