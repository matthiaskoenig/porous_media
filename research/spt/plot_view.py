"""Project native extruded SPT meshes with one shared orthographic camera."""

from dataclasses import dataclass
from functools import cached_property

import meshio
import numpy as np

from research.spt.plot_mesh import boundary_edges, boundary_segments, exterior_edges

GEOMETRY_VIEW = {
    "version": 1,
    "projection": "orthographic",
    "tilt_from_top_degrees": 20.0,
    "azimuth_degrees": -60.0,
    "thickness": "native mesh coordinates; no exaggeration",
    "surface_colors": "native cell values; no lighting or point averaging",
}


def camera_matrix() -> np.ndarray:
    """Foreshorten the XY plane without adding an in-plane rotation or perspective."""
    tilt = np.deg2rad(float(GEOMETRY_VIEW["tilt_from_top_degrees"]))
    azimuth = np.deg2rad(float(GEOMETRY_VIEW["azimuth_degrees"]))
    direction = np.array([np.cos(azimuth), np.sin(azimuth)])
    rotation = np.eye(3)
    rotation[:2, :2] += (np.cos(tilt) - 1) * np.outer(direction, direction)
    rotation[:2, 2] = -np.sin(tilt) * direction
    rotation[2, :2] = np.sin(tilt) * direction
    rotation[2, 2] = np.cos(tilt)
    return rotation


def project(xy: np.ndarray, height: np.ndarray | float) -> np.ndarray:
    """Project positions and their true heights into the common camera plane."""
    coordinates = np.concatenate(
        (xy, np.broadcast_to(height, xy.shape[:-1])[..., None]), axis=-1
    )
    return (coordinates @ camera_matrix().T)[..., :2]


def face_heights(points: np.ndarray, cells: list[meshio.CellBlock]) -> np.ndarray:
    """Retain native lower/upper heights in the planar triangulation's face order."""
    heights = []
    for block in cells:
        vertices = points[block.data]
        if block.type in {"wedge", "hexahedron"}:
            count = 3 if block.type == "wedge" else 4
            lower, upper = vertices[:, :count, 2], vertices[:, count:, 2]
            tolerance = max(float(np.ptp(points, axis=0).max()) * 1e-7, 1e-15)
            if np.any(np.ptp(lower, axis=1) > tolerance) or np.any(
                np.ptp(upper, axis=1) > tolerance
            ):
                raise ValueError("Tilted plotting requires horizontal extrusion caps")
            bounds = np.sort(np.column_stack((lower[:, 0], upper[:, 0])), axis=1)
        elif block.type in {"triangle", "quad"}:
            bounds = np.repeat(vertices[:, :1, 2], 2, axis=1)
        else:
            raise ValueError(f"Unsupported tilted cell type: {block.type}")
        heights.extend([bounds] * (2 if block.type in {"quad", "hexahedron"} else 1))
    return np.concatenate(heights)


@dataclass
class SurfaceProjection:
    """Keep projected caps and visible exterior walls tied to native face indices."""

    faces: np.ndarray
    heights: np.ndarray

    @cached_property
    def caps(self) -> np.ndarray:
        """Display each original cell on its upper cap, without scalar interpolation."""
        return project(self.faces, self.heights[:, 1, None])

    @cached_property
    def edges(self) -> tuple[np.ndarray, np.ndarray]:
        """Join matching seam nodes before finding actual exterior walls and holes."""
        return exterior_edges(self.faces)

    @cached_property
    def walls(self) -> tuple[np.ndarray, np.ndarray]:
        """Show only camera-facing exterior walls, painted from far to near."""
        edges, owners = self.edges
        delta = edges[:, 1] - edges[:, 0]
        normal = np.column_stack((delta[:, 1], -delta[:, 0]))
        away = edges.mean(axis=1) - self.faces[owners].mean(axis=1)
        normal[np.einsum("ij,ij->i", normal, away) < 0] *= -1
        visible = (normal @ camera_matrix()[2, :2] > 0) & (
            self.heights[owners, 1] > self.heights[owners, 0]
        )
        edges, owners = edges[visible], owners[visible]
        order = np.argsort(edges.mean(axis=1) @ camera_matrix()[2, :2])
        edges, owners = edges[order], owners[order]
        lower = project(edges, self.heights[owners, 0, None])
        upper = project(edges, self.heights[owners, 1, None])
        return np.concatenate((lower, upper[:, ::-1]), axis=1), owners

    @cached_property
    def polygons(self) -> list[np.ndarray]:
        """Paint walls before caps; all caps lie above their vertical extrusion."""
        return [*self.walls[0], *self.caps]

    @cached_property
    def face_indices(self) -> np.ndarray:
        """Map every painted polygon back to exactly one original planar face."""
        return np.concatenate((self.walls[1], np.arange(len(self.faces))))

    @cached_property
    def bounds(self) -> np.ndarray:
        """Include the entire visible thickness when framing or cropping a figure."""
        vertices = np.concatenate(self.polygons)
        return np.stack((vertices.min(axis=0), vertices.max(axis=0)))

    def boundary(self, groups: np.ndarray | None = None) -> np.ndarray:
        """Project roof outlines and visible wall silhouettes, never internal cells."""
        edges, owners = (
            self.edges if groups is None else boundary_edges(self.faces, groups)
        )
        roof = project(edges, self.heights[owners, 1, None])
        walls, wall_owners = self.walls
        if not len(walls):
            return roof
        return np.concatenate(
            (
                roof,
                boundary_segments(
                    walls, None if groups is None else groups[wall_owners]
                ),
            )
        )
