"""Define SPT scalar fields, cache selections, and default display settings once."""

from dataclasses import dataclass
from typing import Literal


@dataclass(frozen=True)
class Layer:
    """Describe scalar labels, colors and optional explicit display limits."""

    sid: str
    title: str
    colormap: str = "magma"
    limits: tuple[float, float] | None = None
    factor: float = 1.0
    scale: Literal["linear", "symlog"] = "linear"


LAYERS = (
    Layer("rr_(S_ext)", "Substrate [mM]"),
    Layer("rr_(P_ext)", "Product [mM]", "viridis"),
    Layer("rr_protein", "Protein", "viridis"),
    Layer("rr_(T)", "Toxicant [mM]"),
    Layer("rr_necrosis", "Necrosis (fraction)", "Reds", (0.0, 1.0)),
)
SUMMARY_FIELDS = tuple(layer.sid for layer in LAYERS)
VOLUME = "element_volume_point_TPM"
POSITION = "rr_position"
DIAGNOSTICS = (VOLUME, "rr_Vext", "rr_Vli", "pressure")
ANALYSIS_FIELDS = (
    *SUMMARY_FIELDS,
    VOLUME,
    POSITION,
    *DIAGNOSTICS[1:],
    "rr_zonation_pattern",
    "cell_type",
)
