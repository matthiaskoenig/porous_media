"""Helpers for zonated meshes."""

from pathlib import Path
from typing import Callable, Dict, List, Optional

import meshio
import numpy as np

from porous_media import RESOURCES_DIR
from porous_media.console import console
from porous_media.data.xdmf_tools import AttributeType
from porous_media.mesh.mesh_tools import mesh_to_vtk, mesh_to_xdmf
from porous_media.visualization.image_manipulation import merge_images
from porous_media.visualization.pyvista_visualization import (
    DataLayer,
    DataRangeType,
    VisualizationSettings,
    visualize_data_layers,
)


class ZonationPatterns:
    """Definition of standard zonation patterns.

    This calculates zonated variables based on the position of the mesh nodes.
    """

    @staticmethod
    def position(p: np.ndarray, value: float = 1.0) -> np.ndarray:
        """Pattern based on position information."""
        return value * p

    @staticmethod
    def constant(p: np.ndarray, value: float = 0.5) -> np.ndarray:
        """Constant zonation with value."""
        return value * np.ones_like(p)

    @staticmethod
    def random(
        p: np.ndarray, value_min: float = 0.0, value_max: float = 1.0
    ) -> np.ndarray:
        """Random zonation in [min_value, max_value]."""
        return value_min + (value_max - value_min) * np.random.rand(*p.shape)

    @staticmethod
    def linear_increase(
        p: np.ndarray, value_min: float = 0.0, value_max: float = 1.0
    ) -> np.ndarray:
        """Linear increasing pattern in [min_value, max_value]."""
        return value_min + (value_max - value_min) * p

    @staticmethod
    def linear_decrease(
        p: np.ndarray, value_min: float = 0.0, value_max: float = 1.0
    ) -> np.ndarray:
        """Linear decreasing pattern in [min_value, max_value]."""
        return value_min + (value_max - value_min) * (1.0 - p)

    @staticmethod
    def exp_increase(
        p: np.ndarray, value_min: float = 0.0, value_max: float = 1.0
    ) -> np.ndarray:
        """Exponential increasing pattern in [0, 1]."""
        return value_min + (value_max - value_min) * (np.exp(p) - 1.0) / (
            np.exp(1.0) - 1.0
        )

    @staticmethod
    def exp_decrease(
        p: np.ndarray, value_min: float = 0.0, value_max: float = 1.0
    ) -> np.ndarray:
        """Exponential decreasing pattern in [0, 1]."""
        return 1.0 - ZonationPatterns.exp_increase(p, value_min, value_max)

    @staticmethod
    def only_periportal(
        p: np.ndarray, value_min: float = 0.0, value_max: float = 1.0
    ) -> np.ndarray:
        """Sharp periportal pattern in [0, 1]."""
        data = np.zeros_like(p)
        data[p <= 0.2] = value_max
        data[p > 0.2] = value_min
        return data

    @staticmethod
    def only_pericentral(
        p: np.ndarray, value_min: float = 0.0, value_max: float = 1.0
    ) -> np.ndarray:
        """Sharp pericentral pattern in [0, 1]."""
        data = np.zeros_like(p)
        data[p >= 0.8] = value_max
        data[p < 0.2] = value_min
        return data

    @staticmethod
    def sharp_periportal(
        p: np.ndarray, value_min: float = 0.0, value_max: float = 1.0, n: float = 10.0
    ) -> np.ndarray:
        """Sharp periportal pattern in [0, 1]."""
        return value_min + (value_max - value_min) * (1 - p**n / (p**n + 0.25**n))

    @staticmethod
    def sharp_pericentral(
        p: np.ndarray, value_min: float = 0.0, value_max: float = 1.0, n: float = 10.0
    ) -> np.ndarray:
        """Sharp pericentral pattern in [0, 1]."""
        return value_min + (value_max - value_min) * p**n / (p**n + 0.75**n)


class ZonatedMesh:
    """Class for zonated meshes."""

    patterns: List[Callable] = [
        ZonationPatterns.position,
        ZonationPatterns.constant,
        ZonationPatterns.random,
        ZonationPatterns.linear_increase,
        ZonationPatterns.linear_decrease,
        ZonationPatterns.exp_increase,
        ZonationPatterns.sharp_pericentral,
        ZonationPatterns.only_pericentral,
        ZonationPatterns.exp_decrease,
        ZonationPatterns.sharp_periportal,
        ZonationPatterns.only_periportal,
    ]

    @classmethod
    def create_zonated_mesh_from_vtk(
        cls,
        vtk_path: Path,
        patterns: Optional[List[Callable]] = None,
        remove_point_data: bool = True,
        remove_cell_data: bool = True,
        copy_mesh: bool = True,
    ) -> meshio.Mesh:
        """Create a zonated mesh from VTK and serialize results to XDMF."""
        # patterns
        f_patterns: List[Callable] = cls.patterns if patterns is None else patterns

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
                mesh=m,
                variable_id=f"pattern__{f_pattern.__name__}",
                f_zonation=f_pattern,
            )

        return m

    @staticmethod
    def create_zonated_mesh(
        mesh: meshio.Mesh,
        remove_point_data: bool = True,
        remove_cell_data: bool = True,
        copy_mesh: bool = True,
    ) -> meshio.Mesh:
        """Calculate the distance from periportal and perivenous.

        Uses the cell_type variable for determining the position in the lobulus (periportal, perivenous).
        'cell_type':
            0: internal node
            1: periportal (inflow)
            2: perivenous (outflow)
        """

        # clone mesh
        if copy_mesh:
            m: meshio.Mesh = mesh.copy()
        else:
            m = mesh

        # check for variables
        if "cell_type" not in m.cell_data:
            raise IOError("'cell_type' required in cell_data for zonation patterns")

        # remove cell data
        if remove_cell_data:
            m.cell_data = {"cell_type": m.cell_data["cell_type"]}

        # remove point data
        if remove_point_data:
            m.point_data = {}

        # calculate positions based on cell_type info
        cell_type: np.ndarray = m.cell_data["cell_type"][0]
        position = np.NaN * np.ones_like(cell_type)

        pp_cells: Dict[int, np.ndarray] = {}
        pv_cells: Dict[int, np.ndarray] = {}
        inner_cells: Dict[int, np.ndarray] = {}
        cell_block: meshio._mesh.CellBlock

        for _, cell_block in enumerate(m.cells):
            for kc, cell in enumerate(cell_block.data):
                # calculate the center of mass of cell
                count = 0
                center = np.zeros(shape=3)
                for kp in cell:
                    center += m.points[kp]
                    count += 1
                center = center / count

                # periportal cell
                if np.isclose(cell_type[kc], 1.0):
                    position[kc] = 0
                    pp_cells[kc] = center

                # pericentral cell
                elif np.isclose(cell_type[kc], 2.0):
                    position[kc] = 1.0
                    pv_cells[kc] = center

                # inner cell
                elif np.isclose(cell_type[kc], 0.0):
                    position[kc] = np.NaN
                    inner_cells[kc] = center

        # calculate pp-pv position for all inner cells
        for kc, center in inner_cells.items():
            # shortest distance periportal
            dpv = 1.0
            for center_pv in pv_cells.values():
                d = np.linalg.norm(center - center_pv)
                if d < dpv:
                    dpv = d

            # shortest distance pericentral
            dpp = 1.0
            for center_pp in pp_cells.values():
                d = np.linalg.norm(center - center_pp)
                if d < dpp:
                    dpp = d

            position[kc] = dpp / (dpv + dpp)

        m.cell_data["position"] = [position]
        return m

    @staticmethod
    def _add_zonated_variable(
        mesh: meshio.Mesh, variable_id: str, f_zonation: Callable
    ) -> meshio.Mesh:
        """Add zonation variable to the mesh based on the position variable in [0,1].

        :param variable_id: identifier in cell_data
        :param f_zonation: function to calculate zonation

        position:
            0: periportal
            1: perivenous/pericentral
        """
        # check for variables
        m = mesh
        if "position" not in m.cell_data:
            raise IOError(
                "'position' required in calculate zonation patterns, 'create_zonated_mesh' first."
            )

        position: np.ndarray = m.cell_data["position"][0]
        data = f_zonation(position)
        m.cell_data[variable_id] = [data]
        return m


data_layers = []

for key in [p.__name__ for p in ZonatedMesh.patterns]:
    data_layers.append(
        DataLayer(
            sid=f"pattern__{key}",
            title=key,
            viz_type=AttributeType.SCALAR,
            colormap="magma",
        )
    )

console.print(data_layers)


def visualize_zonation_patterns(
    mesh: meshio.Mesh,
    results_path: Path,
    image_name: str,
    data_layers: List[DataLayer],
    direction: str = "horizontal",
    drange_type: DataRangeType = DataRangeType.LOCAL,
) -> None:
    """Visualize zonation patterns."""
    dl_sids: List[str] = [dl.sid for dl in data_layers]

    # calculate the data ranges
    dmins: Dict[str, float] = {}
    dmaxs: Dict[str, float] = {}

    for key in dl_sids:
        data = mesh.cell_data[key][0]
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
    images: List[Path] = []
    for scalar in data_layers:
        img_path = results_path / scalar.sid / f"{image_name}.png"
        images.append(img_path)

    image: Path = results_path / f"{image_name}.png"
    merge_images(paths=images, direction=direction, output_path=image)
    console.print(f"Image created: file://{image}")


def example_mesh_zonation(results_dir: Path, visualize: bool = True) -> None:
    """Run example for mesh zonation."""

    results_path: Path = results_dir / "mesh_zonation"
    results_path.mkdir(parents=True, exist_ok=True)
    vtk_path = RESOURCES_DIR / "zonation" / "mesh_zonation.vtk"

    # create zonated mesh
    zm = ZonatedMesh()
    m: meshio.Mesh = zm.create_zonated_mesh_from_vtk(vtk_path=vtk_path)

    # serialize mesh with results
    xdmf_path = results_path / "mesh_zonation.xdmf"
    mesh_to_xdmf(m=m, xdmf_path=xdmf_path)
    vtk_path = results_path / "mesh_zonation.vtk"
    mesh_to_vtk(m=m, vtk_path=vtk_path)

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
        "pattern__sharp_pericentral",
        # "pattern__only_pericentral",
        "pattern__linear_decrease",
        "pattern__sharp_periportal",
        # "pattern__only_periportal",
    ]
    data_layers_dict = {dl.sid: dl for dl in data_layers}
    data_layers_selected = [data_layers_dict[key] for key in selected_patterns]
    visualize_zonation_patterns(
        m,
        results_path=results_path,
        image_name="mesh_zonation_selected",
        direction="horizontal",
        data_layers=data_layers_selected,
        drange_type=DataRangeType.GLOBAL,
    )


if __name__ == "__main__":
    from porous_media import RESULTS_DIR

    example_mesh_zonation(results_dir=RESULTS_DIR)
