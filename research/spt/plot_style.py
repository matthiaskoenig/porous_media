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
# A study may add a pattern before this repository names it; keep such a run
# plottable and visibly distinct from the six the manuscript describes.
EXTRA_PATTERN_COLORMAPS = ("YlOrBr", "BuGn", "PuRd", "YlGn")
# Keep colors stable when plotting a subset of the archived resolutions.
CONVERGENCE_ORDER = ("50000", "37250", "25000", "18750", "12500", "09375", "06250")


def simulation_color(pattern: int, dose: int, count: int = len(boundary_flows)) -> str:
    """Return the original pattern palette and one-based dose shade."""
    import matplotlib as mpl
    from matplotlib.colors import to_hex

    name = PATTERN_COLORMAPS.get(
        pattern,
        EXTRA_PATTERN_COLORMAPS[
            (pattern - len(PATTERN_COLORMAPS)) % len(EXTRA_PATTERN_COLORMAPS)
        ],
    )
    return to_hex(mpl.colormaps[name](dose / count), keep_alpha=True)


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
