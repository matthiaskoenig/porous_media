"""Tile native SPT lobule faces into a square without changing their cell values."""

from dataclasses import dataclass

import numpy as np

from research.spt.plot_mesh import boundary_segments


def polygon_area(vertices: np.ndarray) -> float:
    """Return an unsigned polygon area without assuming a vertex orientation."""
    return (
        abs(
            float(
                np.sum(
                    vertices[:, 0] * np.roll(vertices[:, 1], -1)
                    - vertices[:, 1] * np.roll(vertices[:, 0], -1)
                )
            )
        )
        / 2
    )


def clip_square(vertices: np.ndarray, half_width: float) -> np.ndarray:
    """Clip a convex face exactly to four square edges, retaining partial cells."""
    result = vertices
    for dimension, sign in ((0, 1), (0, -1), (1, 1), (1, -1)):
        if not len(result):
            break
        output = []
        previous = result[-1]
        previous_inside = sign * previous[dimension] <= half_width
        for current in result:
            inside = sign * current[dimension] <= half_width
            if inside != previous_inside:
                weight = (sign * half_width - previous[dimension]) / (
                    current[dimension] - previous[dimension]
                )
                crossing = previous + weight * (current - previous)
                crossing[dimension] = sign * half_width
                output.append(crossing)
            if inside:
                output.append(current)
            previous, previous_inside = current, inside
        result = np.asarray(output).reshape(-1, 2)
    return result


def lobule_width(faces: np.ndarray) -> float:
    """Fit the six native outer planes; width is the distance between flat sides.

    SPT lobules are centered at the origin with a pair of vertical flat sides.
    Vessel notches truncate the tips, so their bounding-box height cannot define
    the hexagonal lattice. Flat boundary segments determine the actual spacing.
    """
    if faces.ndim != 3 or faces.shape[1:] != (3, 2) or not np.isfinite(faces).all():
        raise ValueError("Tissue plotting requires finite native XY triangles")
    edges = boundary_segments(faces)
    delta = edges[:, 1] - edges[:, 0]
    lengths = np.linalg.norm(delta, axis=1)
    centers = edges.mean(axis=1)
    supports = []
    for angle in np.arange(6) * np.pi / 3:
        normal = np.array([np.cos(angle), np.sin(angle)])
        distance = centers @ normal
        parallel = np.abs(delta @ normal) <= lengths * 2e-4
        selected = parallel & (distance >= 0.95 * distance.max())
        if not np.any(selected):
            raise ValueError(
                "Tissue plotting requires a complete pointy-top hexagonal lobule"
            )
        supports.append(float(np.median(distance[selected])))
    radius = float(np.median(supports))
    if radius <= 0 or not np.allclose(supports, radius, rtol=1e-4, atol=0):
        raise ValueError(
            "Lobule outer planes must be regular and centered at the origin"
        )
    return 2 * radius


@dataclass(frozen=True)
class TissuePatch:
    """Store clipped native faces, owning cells and lattice translations in lobule widths."""

    faces: np.ndarray
    cell_ids: np.ndarray
    lobule_ids: np.ndarray
    centers: np.ndarray
    width: float
    source_width: float
    nominal_area: float

    @property
    def tissue_area(self) -> float:
        """Measure actual mesh coverage, preserving the source's vascular openings."""
        first, second = (
            self.faces[:, 1] - self.faces[:, 0],
            self.faces[:, 2] - self.faces[:, 0],
        )
        return float(
            np.abs(first[:, 0] * second[:, 1] - first[:, 1] * second[:, 0]).sum() / 2
        )


def tile_lobules(
    faces: np.ndarray, cell_ids: np.ndarray, width: float = 2.5
) -> TissuePatch:
    """Cover a square with repeated lobules and clip every intersecting native cell.

    Coordinates are normalized by the native flat-to-flat lobule width. The
    lattice covers the whole square, including its corners. Native vessel lumens
    remain empty; no scalar values are invented or interpolated between lobules.
    """
    if not np.isfinite(width) or width < 2.5:
        raise ValueError("Square width must be finite and at least 2.5 lobules")
    source_width = lobule_width(faces)
    if (
        cell_ids.shape != (len(faces),)
        or not np.issubdtype(cell_ids.dtype, np.integer)
        or np.any(cell_ids < 0)
    ):
        raise ValueError("Each native face must retain a nonnegative integer cell ID")
    native = faces / source_width
    half = width / 2
    radius = 1 / np.sqrt(3)
    angles = np.pi / 6 + np.arange(6) * np.pi / 3
    hexagon = radius * np.column_stack((np.cos(angles), np.sin(angles)))
    reach = int(np.ceil(width + 2))
    pieces, ids, owners, centers = [], [], [], []
    nominal_area = 0.0
    for row in range(-reach, reach + 1):
        for column in range(-reach, reach + 1):
            center = np.array([column + row / 2, row * np.sqrt(3) / 2])
            hex_clip = clip_square(hexagon + center, half)
            if len(hex_clip) < 3 or polygon_area(hex_clip) < 1e-12:
                continue
            nominal_area += polygon_area(hex_clip)
            lobule = len(centers)
            centers.append(center)
            translated = native + center
            low, high = translated.min(axis=1), translated.max(axis=1)
            intersects = ((low < half) & (high > -half)).all(axis=1)
            inside = ((low >= -half) & (high <= half)).all(axis=1)
            pieces.extend(translated[inside])
            ids.extend(cell_ids[inside])
            owners.extend([lobule] * int(inside.sum()))
            for index in np.flatnonzero(intersects & ~inside):
                clipped = clip_square(translated[index], half)
                if len(clipped) < 3:
                    continue
                for vertex in range(1, len(clipped) - 1):
                    triangle = clipped[[0, vertex, vertex + 1]]
                    if polygon_area(triangle) > 1e-14:
                        pieces.append(triangle)
                        ids.append(cell_ids[index])
                        owners.append(lobule)
    if not np.isclose(nominal_area, width**2, rtol=1e-12, atol=1e-12):
        raise ValueError(
            "Lobule lattice does not completely cover the requested square"
        )
    return TissuePatch(
        np.asarray(pieces),
        np.asarray(ids, dtype=np.int64),
        np.asarray(owners, dtype=np.int64),
        np.asarray(centers),
        width,
        source_width,
        nominal_area,
    )
