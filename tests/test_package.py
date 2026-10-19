"""Verify metadata and resources in editable installs and installed wheels."""

from importlib.metadata import distribution, metadata, version

from porous_media import RESOURCES_DIR, __version__


def test_distribution_metadata() -> None:
    """Expose accurate package metadata to package indexes and installers."""
    assert version("porous-media") == __version__
    project = metadata("porous-media")
    assert "porous media" in project["Summary"]
    assert any("Documentation" in url for url in project.get_all("Project-URL", []))


def test_distributed_resources() -> None:
    """Ship the files needed by the documented examples, including HDF data."""
    for relative_path in (
        "vtk/vtk_single/lobule_zonation_pattern.t0022.vtk",
        "zonation/mesh_zonation.xdmf",
        "zonation/mesh_zonation.h5",
        "zonation/mesh_zonation_lobulus.vtk",
    ):
        assert (RESOURCES_DIR / relative_path).is_file(), relative_path


def test_distribution_excludes_research_code_and_dependencies() -> None:
    """Install the reusable library without publishing project workflows or their tools."""
    package = distribution("porous-media")
    assert all(
        "analyses/" not in str(path) and not str(path).startswith("research/")
        for path in package.files or []
    )
    requirements = package.metadata.get_all("Requires-Dist", [])
    assert not any(
        requirement.lower().startswith(("pyarrow", "openpyxl", "pandas>"))
        for requirement in requirements
    )
