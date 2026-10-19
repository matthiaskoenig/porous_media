"""Definition of SPT information."""

from __future__ import annotations

import os
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from porous_media.visualization.pyvista_visualization import DataLayer
from research import DATA_DIR

data_spt_dir = Path(
    os.environ.get("POROUS_MEDIA_SPT_DATA", DATA_DIR / "spt")
).expanduser()

# subset of scalars to visualize
selection_spt: list[str] = [
    "rr_(S_ext)",
    "rr_(P_ext)",
    "rr_protein",
    "rr_(T)",
    "rr_necrosis",
    # "pressure",
]


def __getattr__(name: str) -> list[DataLayer]:
    """Load legacy visualization exports only when a rendering script requests them."""
    if name == "data_layers_spt":
        from research.spt.visualization import data_layers_spt

        return data_layers_spt
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
