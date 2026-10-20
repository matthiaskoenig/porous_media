"""Extract planar domain outlines without convex hulls or simulation averaging."""

import numpy as np
from scipy.spatial import KDTree

OUTLINE_VERSION = 1


def boundary_segments(
    faces: np.ndarray, groups: np.ndarray | None = None
) -> np.ndarray:
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
    return boundary
