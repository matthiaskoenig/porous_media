"""Extract planar domain outlines without convex hulls or simulation averaging."""

import numpy as np
from scipy.spatial import KDTree

OUTLINE_VERSION = 1


def boundary_segments(
    faces: np.ndarray, groups: np.ndarray | None = None
) -> np.ndarray:
    """Return external edges, preserving concavities, holes and optional sectors."""
    return boundary_edges(faces, groups)[0]


def boundary_edges(
    faces: np.ndarray, groups: np.ndarray | None = None
) -> tuple[np.ndarray, np.ndarray]:
    """Return edges belonging to one face, optionally within independent sectors.

    Coordinate matching joins reflected copies of seam nodes despite roundoff.
    Boundary pairs within the archived reflection precision are also joined,
    with tolerance bounded by 0.1% of the shortest boundary edge.
    Keeping incidence counts preserves concavities and holes. No hull or smoothing
    changes the actual domain boundary.
    """
    span = float(np.ptp(faces.reshape(-1, 2), axis=0).max())
    if span <= 0 or not np.isfinite(faces).all():
        raise ValueError("Boundary extraction requires finite nondegenerate faces")
    edges = np.stack((faces, np.roll(faces, -1, axis=1)), axis=2).reshape(-1, 2, 2)
    coordinates = np.rint((edges - faces.min(axis=(0, 1))) / (span * 1e-10)).astype(
        np.int64
    )
    reverse = (coordinates[:, 0, 0] > coordinates[:, 1, 0]) | (
        (coordinates[:, 0, 0] == coordinates[:, 1, 0])
        & (coordinates[:, 0, 1] > coordinates[:, 1, 1])
    )
    coordinates[reverse] = coordinates[reverse, ::-1]
    keys = coordinates.reshape(-1, 4)
    if groups is not None:
        keys = np.column_stack((keys, np.repeat(groups, faces.shape[1])))
    _, indices, counts = np.unique(keys, axis=0, return_index=True, return_counts=True)
    boundary = edges[indices[counts == 1]]
    owners = indices[counts == 1] // faces.shape[1]
    if groups is None and len(boundary) > 1:
        # Archived sector transforms differ by a few nanometres at their seams.
        # Drop coincident edge pairs without moving vertices or smoothing holes.
        lengths = np.linalg.norm(boundary[:, 1] - boundary[:, 0], axis=1)
        tolerance = min(span * 1e-5, float(lengths.min()) * 1e-3)
        pairs = KDTree(boundary.mean(axis=1)).query_pairs(
            tolerance, output_type="ndarray"
        )
        first, second = boundary[pairs[:, 0]], boundary[pairs[:, 1]]
        coincident = (np.linalg.norm(first - second, axis=2) <= tolerance).all(
            axis=1
        ) | (np.linalg.norm(first - second[:, ::-1], axis=2) <= tolerance).all(axis=1)
        keep = np.ones(len(boundary), dtype=bool)
        keep[pairs[coincident].ravel()] = False
        boundary = boundary[keep]
        owners = owners[keep]
    return boundary, owners


def exterior_edges(faces: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Remove covered interfaces even when neighboring meshes divide edges differently.

    Incidence matching alone leaves T-junctions on mixed-resolution sector seams.
    Remove only overlapping, collinear intervals with opposing outward normals.
    Retained partial edges keep their original owning face for scalar coloring.
    """
    edges, owners = boundary_edges(faces)
    delta = edges[:, 1] - edges[:, 0]
    lengths = np.linalg.norm(delta, axis=1)
    span = float(np.ptp(faces.reshape(-1, 2), axis=0).max())
    tolerance = min(span * 1e-5, float(lengths.min()) * 1e-3)
    unit = delta / lengths[:, None]
    normals = np.column_stack((unit[:, 1], -unit[:, 0]))
    away = edges.mean(axis=1) - faces[owners].mean(axis=1)
    normals[np.einsum("ij,ij->i", normals, away) < 0] *= -1
    pairs = KDTree(edges.mean(axis=1)).query_pairs(
        float(lengths.max()) + tolerance, output_type="ndarray"
    )
    first, second = pairs[:, 0], pairs[:, 1]
    opposite = np.einsum("ij,ij->i", normals[first], normals[second]) < 0
    relative = edges[second] - edges[first, :1]
    distance = np.abs(
        relative[..., 0] * unit[first, None, 1]
        - relative[..., 1] * unit[first, None, 0]
    )
    pairs = pairs[opposite & (distance <= tolerance).all(axis=1)]
    covered: dict[int, list[tuple[float, float]]] = {}
    for first, second in pairs:
        for target, other in ((first, second), (second, first)):
            positions = (edges[other] - edges[target, 0]) @ unit[target]
            low, high = (
                max(0.0, float(positions.min())),
                min(float(lengths[target]), float(positions.max())),
            )
            if high - low > tolerance:
                covered.setdefault(int(target), []).append((low, high))
    if not covered:
        return edges, owners
    retained, retained_owners = [], []
    for index, edge in enumerate(edges):
        start = 0.0
        intervals = [
            *sorted(covered.get(index, [])),
            (float(lengths[index]), float(lengths[index])),
        ]
        for low, high in intervals:
            if low - start > tolerance:
                retained.append(edge[0] + np.array([start, low])[:, None] * unit[index])
                retained_owners.append(owners[index])
            start = max(start, high)
    return np.asarray(retained), np.asarray(retained_owners, dtype=np.int64)
