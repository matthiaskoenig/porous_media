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
    scalar_bar: bool = True
    show_edges: bool = True
    edge_color: str = "darkgray"
    line_width: float = 0.48  # Matplotlib points: one pixel at the default 150 DPI.
    edge_alpha: float = 0.18
    outline: bool = False


SPT_LAYERS = (
    Layer("rr_necrosis", "Necrosis (0: alive, 1: death)", "binary"),
    Layer("rr_protein", "Protein", "hot"),
    Layer("rr_(S_ext)", "Substrate S plasma [mM]"),
    Layer("rr_(P_ext)", "Product P plasma [mM]"),
    Layer("rr_(S)", "Substrate S [mM]"),
    Layer("rr_(P)", "Product P [mM]"),
    Layer("rr_(T)", "Toxic compound T [mM]"),
    Layer("pressure", "Pressure [?]"),
)
SUMMARY_FIELDS = ("rr_(S_ext)", "rr_(P_ext)", "rr_protein", "rr_(T)", "rr_necrosis")
LAYERS = tuple(
    next(layer for layer in SPT_LAYERS if layer.sid == sid) for sid in SUMMARY_FIELDS
)
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
