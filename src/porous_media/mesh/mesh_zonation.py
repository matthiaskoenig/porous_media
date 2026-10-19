"""Helpers for zonated meshes."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, ClassVar, Protocol

import meshio
import numpy as np
from scipy.spatial import KDTree

from porous_media import RESOURCES_DIR
from porous_media.console import console
from porous_media.data.xdmf_tools import AttributeType
from porous_media.mesh.mesh_tools import mesh_to_vtk, mesh_to_xdmf

if TYPE_CHECKING:
    from porous_media.visualization.pyvista_visualization import (
        DataLayer,
        DataRangeType,
    )


class ZonationPattern(Protocol):
    """Named function that creates zonation values and accepts a scaling factor."""

    @property
    def __name__(self) -> str:
        """Return the name used to identify the resulting variable."""
        ...

    def __call__(self, p: np.ndarray, *, f_scale: float = 1.0) -> np.ndarray:
        """Calculate the pattern at the supplied positions."""
        ...


class ZonationPatterns:
    """Definition of standard zonation patterns.

    This calculates zonated variables based on the position of the mesh nodes.
    """

    @staticmethod
    def position(p: np.ndarray, value: float = 1.0, f_scale: float = 1.0) -> np.ndarray:
        """Pattern based on position information."""
        return f_scale * value * p

    @staticmethod
    def constant(p: np.ndarray, value: float = 0.5, f_scale: float = 1.0) -> np.ndarray:
        """Constant zonation with value."""
        return f_scale * value * np.ones_like(p)

    @staticmethod
    def random(
        p: np.ndarray,
        value_min: float = 0.0,
        value_max: float = 1.0,
        f_scale: float = 1.0,
        *,
        rng: np.random.Generator | None = None,
    ) -> np.ndarray:
        """Draw random zonation using a local or caller-supplied generator."""
        generator = np.random.default_rng() if rng is None else rng
        return f_scale * (
            value_min + (value_max - value_min) * generator.random(p.shape)
        )

    @staticmethod
    def linear_increase(
        p: np.ndarray,
        value_min: float = 0.0,
        value_max: float = 1.0,
        f_scale: float = 1.0,
    ) -> np.ndarray:
        """Linear increasing pattern in [min_value, max_value]."""
        return f_scale * (value_min + (value_max - value_min) * p)

    @staticmethod
    def linear_decrease(
        p: np.ndarray,
        value_min: float = 0.0,
        value_max: float = 1.0,
        f_scale: float = 1.0,
    ) -> np.ndarray:
        """Linear decreasing pattern in [min_value, max_value]."""
        return f_scale * (value_min + (value_max - value_min) * (1.0 - p))

    @staticmethod
    def exp_increase(
        p: np.ndarray,
        value_min: float = 0.0,
        value_max: float = 1.0,
        f_scale: float = 1.0,
    ) -> np.ndarray:
        """Exponential increasing pattern in [0, 1]."""
        result: np.ndarray = f_scale * (
            value_min
            + (value_max - value_min) * (np.exp(p) - 1.0) / (np.exp(1.0) - 1.0)
        )
        return result

    @staticmethod
    def exp_decrease(
        p: np.ndarray,
        value_min: float = 0.0,
        value_max: float = 1.0,
        f_scale: float = 1.0,
    ) -> np.ndarray:
        """Exponential decreasing pattern in [0, 1]."""
        return f_scale * (
            value_min
            + value_max
            - ZonationPatterns.exp_increase(p, value_min, value_max)
        )

    @staticmethod
    def only_periveneous(
        p: np.ndarray,
        value_min: float = 0.0,
        value_max: float = 1.0,
        f_scale: float = 1.0,
    ) -> np.ndarray:
        """Only periveneous pattern in [0, 1]."""
        data = np.full_like(p, value_min, dtype=float)
        data[p >= 0.8] = value_max
        data[p < 0.2] = value_min
        return f_scale * data

    @staticmethod
    def only_periportal(
        p: np.ndarray,
        value_min: float = 0.0,
        value_max: float = 1.0,
        f_scale: float = 1.0,
    ) -> np.ndarray:
        """Only periportal pattern in [0, 1]."""
        data = np.full_like(p, value_min, dtype=float)
        data[p <= 0.2] = value_max
        data[p > 0.2] = value_min
        return f_scale * data

    @staticmethod
    def sharp_periportal(
        p: np.ndarray,
        value_min: float = 0.0,
        value_max: float = 1.0,
        n: float = 10.0,
        f_scale: float = 1.0,
    ) -> np.ndarray:
        """Sharp periportal pattern in [0, 1]."""
        result: np.ndarray = f_scale * (
            value_min + (value_max - value_min) * (1 - p**n / (p**n + 0.25**n))
        )
        return result

    @staticmethod
    def sharp_periveneous(
        p: np.ndarray,
        value_min: float = 0.0,
        value_max: float = 1.0,
        n: float = 10.0,
        f_scale: float = 1.0,
    ) -> np.ndarray:
        """Sharp periveneous pattern in [0, 1]."""
        result: np.ndarray = f_scale * (
            value_min
            + (value_max - value_min) * (1 - (1 - p) ** n / ((1 - p) ** n + 0.25**n))
        )
        return result


class ZonatedMesh:
    """Class for zonated meshes."""

    patterns: ClassVar[list[ZonationPattern]] = [
        # increasing patterns
        ZonationPatterns.linear_increase,
        ZonationPatterns.exp_increase,
        ZonationPatterns.sharp_periveneous,
        ZonationPatterns.only_periveneous,
        # decreasing patterns
        ZonationPatterns.linear_decrease,
        ZonationPatterns.exp_decrease,
        ZonationPatterns.sharp_periportal,
        ZonationPatterns.only_periportal,
        # misc patterns
        ZonationPatterns.position,
        ZonationPatterns.constant,
        ZonationPatterns.random,
    ]

    @classmethod
    def create_zonated_mesh_from_vtk(
        cls,
        vtk_path: Path,
        patterns: list[ZonationPattern] | None = None,
        remove_point_data: bool = True,
        remove_cell_data: bool = True,
        copy_mesh: bool = True,
        *,
        rng: np.random.Generator | None = None,
    ) -> meshio.Mesh:
        """Return a zonated mesh; use ``rng`` for reproducible random patterns."""
        f_patterns = cls.patterns if patterns is None else patterns

        # mesh
        mesh: meshio.Mesh = meshio.read(vtk_path)
        m: meshio.Mesh = cls.create_zonated_mesh(
            mesh,
            remove_point_data=remove_point_data,
            remove_cell_data=remove_cell_data,
            copy_mesh=copy_mesh,
        )

        for f_pattern in f_patterns:
            cls._add_zonated_variable(
                m=m,
                variable_id=f"pattern__{f_pattern.__name__}",
                f_zonation=f_pattern,
                rng=rng,
            )

        return m

    @staticmethod
    def create_zonated_mesh(
        mesh: meshio.Mesh,
        remove_point_data: bool = True,
        remove_cell_data: bool = True,
        copy_mesh: bool = True,
    ) -> meshio.Mesh:
        """Calculate positions using nearest boundary centers across all cell blocks.

        ``cell_type`` labels interior, periportal, and perivenous cells as 0, 1,
        and 2. Both boundary types are required. Supply finite nonnegative volumes
        as ``volume`` or ``element_volume_point_TPM``, with a positive total.
        Geometry and fields are validated before modifying a mesh in place.
        """
        volume_key = (
            "volume" if "volume" in mesh.cell_data else "element_volume_point_TPM"
        )
        for name in ("cell_type", volume_key):
            if name not in mesh.cell_data:
                raise ValueError(f"Cell field {name!r} is required for zonation")
            if len(mesh.cell_data[name]) != len(mesh.cells):
                raise ValueError(f"Cell field {name!r} must match every cell block")
            for block, values in zip(mesh.cells, mesh.cell_data[name], strict=True):
                if np.asarray(values).size != len(block):
                    raise ValueError(
                        f"Cell field {name!r} must have one value per cell"
                    )
        if not mesh.cells or not any(len(block) for block in mesh.cells):
            raise ValueError("Zonation requires at least one cell")
        points = np.asarray(mesh.points, dtype=float)
        if not np.isfinite(points).all():
            raise ValueError("Mesh points must be finite")
        centers = np.concatenate(
            [points[block.data].mean(axis=1) for block in mesh.cells], axis=0
        )
        labels = np.concatenate(mesh.cell_data["cell_type"]).reshape(-1)
        volume = np.concatenate(mesh.cell_data[volume_key]).reshape(-1)
        if not np.isfinite(volume).all() or (volume < 0).any() or volume.sum() <= 0:
            raise ValueError(
                "Cell volumes must be finite, nonnegative, and have a positive total"
            )
        inner, periportal, perivenous = [
            np.isclose(labels, value) for value in (0, 1, 2)
        ]
        if not (inner | periportal | perivenous).all():
            raise ValueError("cell_type values must be 0, 1, or 2")
        if not periportal.any() or not perivenous.any():
            raise ValueError(
                "Both periportal and perivenous boundary cells are required"
            )
        position = np.empty(len(centers))
        position[periportal] = 0.0
        position[perivenous] = 1.0
        if inner.any():
            dpp = KDTree(centers[periportal]).query(centers[inner])[0]
            dpv = KDTree(centers[perivenous]).query(centers[inner])[0]
            distance = dpp + dpv
            if not np.isfinite(distance).all() or (distance <= 0).any():
                raise ValueError(
                    "Interior cells need a finite, nonzero boundary distance"
                )
            position[inner] = dpp / distance
        m = mesh.copy() if copy_mesh else mesh
        if remove_point_data:
            m.point_data = {}
        if remove_cell_data:
            m.cell_data = {
                "cell_type": m.cell_data["cell_type"],
                "volume": m.cell_data[volume_key],
            }
        else:
            m.cell_data["volume"] = m.cell_data[volume_key]
        m.cell_data["position"] = list(
            np.split(position, np.cumsum([len(block) for block in m.cells])[:-1])
        )
        return m

    @staticmethod
    def _add_zonated_variable(
        m: meshio.Mesh,
        variable_id: str,
        f_zonation: ZonationPattern,
        *,
        rng: np.random.Generator | None = None,
    ) -> meshio.Mesh:
        """Sample a pattern once and normalize its volume-weighted mean to 0.5."""
        position = np.concatenate(m.cell_data["position"]).reshape(-1)
        volume = np.concatenate(m.cell_data["volume"]).reshape(-1)
        protein = np.asarray(
            ZonationPatterns.random(position, rng=rng)
            if f_zonation is ZonationPatterns.random
            else f_zonation(position),
            dtype=float,
        )
        if protein.shape != position.shape or not np.isfinite(protein).all():
            raise ValueError("A zonation pattern must return one finite value per cell")
        total = float(np.dot(protein, volume))
        if not np.isfinite(total) or total == 0:
            raise ValueError(
                "A zonation pattern must have a finite, nonzero weighted total"
            )
        scaled = protein * (0.5 * volume.sum() / total)
        m.cell_data[variable_id] = list(
            np.split(scaled, np.cumsum([len(block) for block in m.cells])[:-1])
        )
        return m


def visualize_zonation_patterns(
    mesh: meshio.Mesh,
    results_path: Path,
    image_name: str,
    data_layers: list[DataLayer],
    direction: str = "horizontal",
    drange_type: DataRangeType | None = None,
) -> None:
    """Visualize patterns with local ranges by default, importing rendering lazily."""
    from porous_media.visualization.image_manipulation import merge_images
    from porous_media.visualization.pyvista_visualization import (
        DataRangeType,
        VisualizationSettings,
        visualize_data_layers,
    )

    dl_sids: list[str] = [dl.sid for dl in data_layers]

    # calculate the data ranges
    dmins: dict[str, float] = {}
    dmaxs: dict[str, float] = {}

    for key in dl_sids:
        data = np.concatenate(mesh.cell_data[key])
        # min, max
        dmins[key] = data.min()
        dmaxs[key] = data.max()

    # handle the case of global data ranges
    if drange_type == DataRangeType.GLOBAL:
        dmin_global = min(dmins.values())
        dmax_global = max(dmaxs.values())
        for key in mesh.cell_data:
            dmins[key] = dmin_global
            dmaxs[key] = dmax_global

    for dl in data_layers:
        sid = dl.sid
        dl.color_limits = (dmins[sid], dmaxs[sid])

    # create raw images
    visualize_data_layers(
        mesh=mesh,
        data_layers=data_layers,
        output_dir=results_path,
        image_name=image_name,
        visualization_settings=VisualizationSettings(off_screen=True),
    )
    # combine images
    images: list[Path] = []
    for scalar in data_layers:
        img_path = results_path / scalar.sid / f"{image_name}.png"
        images.append(img_path)

    image: Path = results_path / f"{image_name}.png"
    merge_images(paths=images, direction=direction, output_path=image)
    console.print(f"Image created: file://{image}")


def example_mesh_zonation(results_dir: Path, visualize: bool = True) -> None:
    """Export the zonation example, optionally rendering its patterns."""
    results_path: Path = results_dir / "mesh_zonation"
    results_path.mkdir(parents=True, exist_ok=True)
    vtk_path = RESOURCES_DIR / "zonation" / "mesh_zonation_lobulus.vtk"

    # create zonated mesh
    zm = ZonatedMesh()
    m: meshio.Mesh = zm.create_zonated_mesh_from_vtk(vtk_path=vtk_path)

    # serialize mesh with results
    xdmf_path = results_path / "mesh_zonation_lobulus.xdmf"
    mesh_to_xdmf(m=m, xdmf_path=xdmf_path)
    vtk_path = results_path / "mesh_zonation_lobulus.vtk"
    mesh_to_vtk(m=m, vtk_path=vtk_path)

    if not visualize:
        return

    from porous_media.visualization.pyvista_visualization import (
        DataLayer,
        DataRangeType,
    )

    data_layers = [
        DataLayer(
            sid=f"pattern__{pattern.__name__}",
            title=pattern.__name__,
            viz_type=AttributeType.SCALAR,
        )
        for pattern in ZonatedMesh.patterns
    ]

    # visualize mesh
    console.rule(title="Mesh Visualization", style="white")
    visualize_zonation_patterns(
        m,
        results_path=results_path,
        image_name="mesh_zonation_all",
        data_layers=data_layers,
        direction="square",
    )

    # Second mesh
    console.rule(title="Mesh Visualization", style="white")
    selected_patterns = [
        "pattern__constant",
        "pattern__linear_increase",
        "pattern__sharp_periveneous",
        "pattern__linear_decrease",
        "pattern__sharp_periportal",
        "pattern__random",
    ]
    data_layers_dict = {dl.sid: dl for dl in data_layers}
    data_layers_selected = [data_layers_dict[key] for key in selected_patterns]

    for direction in ["horizontal", "vertical"]:
        visualize_zonation_patterns(
            m,
            results_path=results_path,
            image_name=f"mesh_zonation_selected_{direction}",
            direction=direction,
            data_layers=data_layers_selected,
            # drange_type=DataRangeType.GLOBAL,
            drange_type=DataRangeType.LOCAL,
        )
