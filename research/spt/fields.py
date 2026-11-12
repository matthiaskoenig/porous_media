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
    show_edges: bool = False
    edge_color: str = "darkgray"
    line_width: float = 0.48  # Matplotlib points: one pixel at the default 150 DPI.
    edge_alpha: float = 0.18
    outline: bool = True


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

# Cell fields a reconstructed mesh sidecar keeps: every display layer, every
# analysis field, and the categorical fields whose temporal rule lives in
# research.spt.sampling. A source run carries 30 cell fields plus mechanical
# vectors and 3x3 tensors that nothing downstream reads; measured on one 1.76 GB
# sidecar, `stress`, `Lagrange_strain` and `solid_stress_TPM` alone were 45% of
# its storage and the vectors another 16%. Reconstruct with every source field
# through the `--all-fields` option when the full mechanical state is wanted.
GEOMETRY_CELL_FIELDS: tuple[str, ...] = tuple(
    dict.fromkeys(
        (
            *(layer.sid for layer in SPT_LAYERS),
            *ANALYSIS_FIELDS,
            "gauss_points_TPM",
        )
    )
)
# No point field is read back; the composition writes its own ``sector_id``.
GEOMETRY_POINT_FIELDS: tuple[str, ...] = ()
