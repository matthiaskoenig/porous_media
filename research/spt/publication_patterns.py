"""Draw reusable monochrome profile icons for manuscript zonation patterns."""

from __future__ import annotations

import numpy as np
from matplotlib.artist import Artist
from matplotlib.backend_bases import RendererBase
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.text import Text
from matplotlib.transforms import Affine2D, Bbox

PATTERN_ICON_DESCRIPTION = (
    "Black-and-white icons show schematic one-dimensional profiles of relative "
    "initial protein abundance (vertical) against periportal-to-perivenous "
    "position (left to right). The random profile is illustrative."
)


def pattern_profile(pattern: int) -> tuple[np.ndarray, np.ndarray]:
    """Return normalized profile paths from periportal to perivenous position.

    Repeated positions retain sharp steps. Profiles explain pattern identities,
    rather than quantitative values or a simulation's particular random cells.
    """
    match pattern:
        case 0:
            position, values = [0, 1], [0.5, 0.5]
        case 1:
            position, values = [0, 1], [0, 1]
        case 2:
            position, values = [0, 0.7, 0.7, 1], [0, 0, 1, 1]
        case 3:
            position, values = [0, 1], [1, 0]
        case 4:
            position, values = [0, 0.3, 0.3, 1], [1, 1, 0, 0]
        case 5:
            edges = np.linspace(0, 1, 10)
            abundance = np.random.default_rng(1729).uniform(0.1, 0.9, 9)
            return np.repeat(edges, 2)[1:-1], np.repeat(abundance, 2)
        case _:
            raise ValueError(f"Unknown zonation pattern: {pattern}")
    return np.asarray(position, dtype=float), np.asarray(values, dtype=float)


class PatternIcon(Artist):
    """Place an upright black-and-white profile panel beside a pattern label."""

    def __init__(self, target: Text, pattern: int, *, size: float = 22) -> None:
        """Anchor to the label's left edge; panel height and gap are in points."""
        super().__init__()
        self.target = target
        self.size = size
        self.gap = 5.0
        position, values = pattern_profile(pattern)
        self.panel = Rectangle(
            (0, 0), 1, 1, facecolor="white", edgecolor="black", linewidth=0.5
        )
        self.profile = Line2D(
            0.12 + 0.76 * position,
            0.14 + 0.72 * values,
            color="black",
            linewidth=1.2,
            solid_capstyle="butt",
        )
        self.set_gid(f"pattern-icon-{chr(97 + pattern)}")
        self.set_in_layout(False)
        self.set_clip_on(False)

    def get_window_extent(self, renderer: RendererBase | None = None) -> Bbox:
        """Include the full profile panel when exports are cropped to content."""
        if renderer is None:
            raise ValueError("Pattern icon bounds require a renderer")
        text = self.target.get_window_extent(renderer)
        height = float(np.asarray(renderer.points_to_pixels(self.size)))
        width = 1.5 * height
        gap = float(np.asarray(renderer.points_to_pixels(self.gap)))
        return Bbox.from_bounds(
            text.x0 - gap - width,
            text.y0 + (text.height - height) / 2,
            width,
            height,
        )

    def draw(self, renderer: RendererBase) -> None:
        """Keep the panel border and profile editable as SVG and PDF vector paths."""
        if not self.get_visible():
            return
        bounds = self.get_window_extent(renderer)
        transform = (
            Affine2D()
            .scale(bounds.width, bounds.height)
            .translate(bounds.x0, bounds.y0)
        )
        parent = self.get_figure(root=True)
        if parent is None:
            raise ValueError("Pattern icon must belong to a figure")
        renderer.open_group("pattern_icon", gid=self.get_gid())
        for artist in (self.panel, self.profile):
            artist.set_figure(parent)
            artist.set_transform(transform)
            artist.set_clip_on(False)
            artist.draw(renderer)
        renderer.close_group("pattern_icon")
        self.stale = False

    def scale(self, factor: float) -> None:
        """Resize the panel, line and gap together with compact-export typography."""
        self.size *= factor
        self.gap *= factor
        for artist in (self.panel, self.profile):
            artist.set_linewidth(artist.get_linewidth() * factor)


def add_pattern_icon(target: Text, pattern: int, *, size: float = 22) -> PatternIcon:
    """Attach the same profile to a heading or a row label in every figure."""
    icon = PatternIcon(target, pattern, size=size)
    if target.axes is not None:
        target.axes.add_artist(icon)
    else:
        parent = target.get_figure()
        if parent is None:
            raise ValueError("Pattern label must belong to a figure")
        parent.add_artist(icon)
    return icon
