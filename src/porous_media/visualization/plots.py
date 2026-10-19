"""Opt-in Matplotlib styling that restores the application's configuration."""

from collections.abc import Generator
from contextlib import contextmanager

import matplotlib

SMALL_SIZE = 13
MEDIUM_SIZE = 18
BIGGER_SIZE = 25


@contextmanager
def plot_style() -> Generator[None]:
    """Apply plot sizes within an explicit context and restore them afterward."""
    with matplotlib.rc_context(
        {
            "font.size": SMALL_SIZE,
            "axes.titlesize": SMALL_SIZE,
            "axes.labelsize": MEDIUM_SIZE,
            "xtick.labelsize": SMALL_SIZE,
            "ytick.labelsize": SMALL_SIZE,
            "legend.fontsize": SMALL_SIZE,
            "figure.titlesize": BIGGER_SIZE,
        }
    ):
        yield
