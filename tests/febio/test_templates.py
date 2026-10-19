"""Verify safe FEBio rendering without external solvers."""

from pathlib import Path
from xml.etree import ElementTree

import pytest
from jinja2 import UndefinedError
from jinja2.sandbox import SecurityError

from porous_media.febio.febio_manipulation import create_feb_files


@pytest.fixture
def template(tmp_path: Path) -> Path:
    """Provide a small XML template with one required parameter."""
    path = tmp_path / "model.feb.template"
    path.write_text('<model value="{{ value }}"/>', encoding="utf-8")
    return path


def test_render_escapes_xml_and_creates_parent(template: Path, tmp_path: Path) -> None:
    """Preserve literal parameter text rather than injecting XML markup."""
    output = tmp_path / "nested" / "inputs"
    value = 'König <mesh> & "parameter"'
    create_feb_files(template, {"simulation": {"value": value}}, output)
    root = ElementTree.parse(output / "model_simulation.feb").getroot()
    assert root.attrib["value"] == value


@pytest.mark.parametrize(
    "sim_id", ["x/../../escaped", "../escaped", "x\\escaped", "", "..", "CON", "a:b"]
)
def test_render_rejects_unsafe_identifiers(
    template: Path, tmp_path: Path, sim_id: str
) -> None:
    """Reject traversal before creating any simulation outputs."""
    output = tmp_path / "outputs"
    (output / "model_x").mkdir(parents=True)
    with pytest.raises(ValueError, match="filename component"):
        create_feb_files(template, {sim_id: {"value": 1}}, output)
    assert not (tmp_path / "escaped.feb").exists()
    assert not list(output.rglob("*.feb"))


def test_render_requires_explicit_overwrite(template: Path, tmp_path: Path) -> None:
    """Preserve an existing model unless replacement is explicitly requested."""
    output = tmp_path / "outputs"
    create_feb_files(template, {"one": {"value": 1}}, output)
    destination = output / "model_one.feb"
    before = destination.read_bytes()
    with pytest.raises(FileExistsError):
        create_feb_files(template, {"one": {"value": 2}}, output)
    assert destination.read_bytes() == before
    create_feb_files(template, {"one": {"value": 2}}, output, overwrite=True)
    assert ElementTree.parse(destination).getroot().attrib["value"] == "2"


def test_render_rejects_python_internal_access(template: Path, tmp_path: Path) -> None:
    """Prevent a template from reaching Python globals through Jinja helpers."""
    template.write_text(
        "{{ cycler.__init__.__globals__.os.getcwd() }}", encoding="utf-8"
    )
    output = tmp_path / "outputs"
    with pytest.raises(SecurityError):
        create_feb_files(template, {"one": {}}, output)
    assert not output.exists()


def test_render_rejects_missing_parameters(template: Path, tmp_path: Path) -> None:
    """Fail clearly instead of silently creating an empty required parameter."""
    output = tmp_path / "outputs"
    with pytest.raises(UndefinedError, match="value"):
        create_feb_files(template, {"one": {}}, output)
    assert not output.exists()


def test_render_preserves_symlink_targets(template: Path, tmp_path: Path) -> None:
    """Never overwrite an unrelated file through a generated output symlink."""
    target = tmp_path / "unrelated.feb"
    target.write_text("unrelated", encoding="utf-8")
    output = tmp_path / "outputs"
    output.mkdir()
    try:
        (output / "model_one.feb").symlink_to(target)
    except OSError:
        pytest.skip("Creating symbolic links is unavailable")
    with pytest.raises(ValueError, match="symbolic link"):
        create_feb_files(template, {"one": {"value": 1}}, output, overwrite=True)
    assert target.read_text(encoding="utf-8") == "unrelated"
