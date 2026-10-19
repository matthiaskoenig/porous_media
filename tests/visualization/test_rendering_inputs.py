"""Check rendering inputs and image workflows without requiring a graphics backend."""

from pathlib import Path
from unittest.mock import MagicMock

import meshio
import numpy as np
import pytest
from PIL import Image

from porous_media.data.xdmf_tools import DataLimits
from porous_media.visualization import pyvista_visualization
from porous_media.visualization.image_annotation import annotate_image_text
from porous_media.visualization.pyvista_visualization import (
    DataLayer,
    create_combined_images,
    visualize_data_layers,
)


@pytest.mark.parametrize(
    ("sid", "image_name"), [("../escape", "image"), ("pressure", "../escape")]
)
def test_rendering_rejects_escaping_paths(
    tmp_path: Path, sid: str, image_name: str
) -> None:
    """Reject path traversal before a plotter or output directory is created."""
    mesh = meshio.Mesh(np.eye(3), [("triangle", [[0, 1, 2]])])
    with pytest.raises(ValueError, match="filename component"):
        visualize_data_layers(
            mesh, [DataLayer(sid, "Pressure")], tmp_path / "output", image_name
        )
    assert not (tmp_path / "output").exists()


def test_combined_images_reuse_generator_selection(tmp_path: Path) -> None:
    """Generate every frame when selection is an iterator."""
    panels = tmp_path / "panels" / "pressure"
    panels.mkdir(parents=True)
    for step in range(2):
        with Image.new("RGB", (10, 10), "red") as image:
            image.save(panels / f"sim_{step:05d}.png")
    outputs = create_combined_images(2, tmp_path, iter(["pressure"]), "horizontal")
    assert len(outputs) == 2
    assert all(path.is_file() for path in outputs)


def test_annotation_uses_portable_font(tmp_path: Path) -> None:
    """Annotate without an operating-system Arial installation and create parents."""
    source = tmp_path / "source.png"
    with Image.new("RGB", (100, 100), "white") as image:
        image.save(source)
    destination = tmp_path / "nested" / "annotation.png"
    annotate_image_text(source, destination, "Text", font_size=20)
    with Image.open(destination) as result:
        assert result.size == (100, 100)
        assert np.asarray(result).min() < 255


def test_rendering_closes_plotter_on_failure(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Release graphics resources when setting up a mesh fails."""
    plotter = MagicMock()
    plotter.add_mesh.side_effect = RuntimeError("render failed")
    monkeypatch.setattr(
        pyvista_visualization.pv, "Plotter", MagicMock(return_value=plotter)
    )
    mesh = meshio.Mesh(
        np.eye(3), [("triangle", [[0, 1, 2]])], cell_data={"pressure": [np.ones(1)]}
    )
    with pytest.raises(RuntimeError, match="render failed"):
        visualize_data_layers(
            mesh, [DataLayer("pressure", "Pressure")], tmp_path, "image"
        )
    plotter.close.assert_called_once()


def test_rendering_converts_values_without_mutating_input(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Apply display units to the rendered mesh while preserving original simulation data."""
    plotter = MagicMock()
    monkeypatch.setattr(
        pyvista_visualization.pv, "Plotter", MagicMock(return_value=plotter)
    )
    mesh = meshio.Mesh(
        np.eye(3), [("triangle", [[0, 1, 2]])], cell_data={"pressure": [np.ones(1)]}
    )
    layer = DataLayer("pressure", "Pressure", conversion_factor=100.0, scalar_bar=False)
    visualize_data_layers(mesh, [layer], tmp_path, "image")
    rendered = plotter.add_mesh.call_args.args[0]
    np.testing.assert_array_equal(rendered.cell_data["pressure"], [100.0])
    np.testing.assert_array_equal(mesh.cell_data["pressure"][0], [1.0])


def test_color_limits_use_display_units() -> None:
    """Keep converted ranges ordered even for a negative conversion factor."""
    layer = DataLayer("pressure", "Pressure", conversion_factor=-2.0)
    layer.update_color_limits(DataLimits({"pressure": (1.0, 3.0)}))
    assert layer.color_limits == (-6.0, -2.0)


def test_vtu_conversion_writes_readable_output(tmp_path: Path) -> None:
    """Create the XDMF file promised by the conversion function's name."""
    from porous_media.visualization.vtu_reader import vtu_to_xdmf

    source = tmp_path / "source.vtu"
    meshio.write(
        source,
        meshio.Mesh(
            np.eye(3),
            [("triangle", [[0, 1, 2]])],
            cell_data={"pressure": [np.array([7.0])]},
        ),
    )
    output = tmp_path / "nested" / "converted.xdmf"
    vtu_to_xdmf(source, output)
    assert output.with_suffix(".h5").is_file()
    result = meshio.read(output)
    np.testing.assert_array_equal(result.cell_data["pressure"][0], [7.0])
    with pytest.raises(FileExistsError):
        vtu_to_xdmf(source, output)
