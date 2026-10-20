"""Preserve the SPT manuscript's plotting conventions across workflows."""

from research.spt.conditions import boundary_flows

LABEL_STYLE = {"fontsize": 12, "fontweight": "bold"}
PATTERN_COLORMAPS = {
    0: "Greys",
    1: "Blues",
    2: "Purples",
    3: "Reds",
    4: "Oranges",
    5: "Greens",
}
# Keep colors stable when plotting a subset of the archived resolutions.
CONVERGENCE_ORDER = ("50000", "37250", "25000", "18750", "12500", "09375", "06250")


def simulation_color(pattern: int, dose: int, count: int = len(boundary_flows)) -> str:
    """Return the original pattern palette and one-based dose shade."""
    import matplotlib as mpl
    from matplotlib.colors import to_hex

    return to_hex(
        mpl.colormaps[PATTERN_COLORMAPS[pattern]](dose / count), keep_alpha=True
    )


def convergence_color(mesh: str, index: int, count: int) -> str:
    """Map coarse to fine meshes onto the manuscript's viridis progression."""
    import matplotlib as mpl
    from matplotlib.colors import to_hex

    fraction = (
        CONVERGENCE_ORDER.index(mesh) / (len(CONVERGENCE_ORDER) - 1)
        if mesh in CONVERGENCE_ORDER
        else index / max(count - 1, 1)
    )
    return to_hex(mpl.colormaps["viridis"](fraction))
