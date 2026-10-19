"""Share temporal interpolation rules across SPT reconstruction and analysis."""

import numpy as np

DISCRETE_FIELDS = frozenset(
    {"cell_type", "rr_zonation_pattern", "rr_necrosis", "gauss_points_TPM"}
)


def bracket(times: np.ndarray, time: float) -> tuple[int, int, float]:
    """Find neighboring native frames without allowing temporal extrapolation."""
    if not np.isfinite(time) or time < times[0] or time > times[-1]:
        raise ValueError("Requested time outside source history")
    upper = int(np.searchsorted(times, time))
    lower = upper if times[upper] == time else upper - 1
    weight = (
        float((time - times[lower]) / (times[upper] - times[lower]))
        if lower != upper
        else 0.0
    )
    return lower, upper, weight


def blend(name: str, low: np.ndarray, high: np.ndarray, weight: float) -> np.ndarray:
    """Interpolate continuous arrays and retain previous discrete/integer values."""
    if low.shape != high.shape:
        raise ValueError(f"Field shape changes: {name}")
    if weight == 0 or name in DISCRETE_FIELDS or low.dtype.kind in "biu":
        return low
    return low * (1.0 - weight) + high * weight
