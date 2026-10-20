"""Protect fast, cached, relocatable report previews and original figure links."""

import shutil
from pathlib import Path

import pytest
from PIL import Image

from research.spt.report import build_report
from research.spt.report_previews import PREVIEW_SIZE, preview_image
from research.tests.test_spt_report import ReportParser, _publish


def test_report_loads_small_previews_and_inspects_originals(tmp_path: Path) -> None:
    """Only small WebP files appear in gallery images; full exports stay on demand."""
    generation = _publish(
        tmp_path, "figures/convergence_sixth/convergence_06250", ("timecourse.png",)
    )
    source = generation / "timecourse.png"
    with Image.new("RGB", (3072, 2048), "navy") as image:
        image.save(source)
    # Republish the output fingerprint after constructing the synthetic large figure.
    import json

    pointer = generation.parent / "manifest.json"
    manifest = json.loads(pointer.read_text())
    manifest["outputs"][0].update(
        size=source.stat().st_size, mtime_ns=source.stat().st_mtime_ns
    )
    pointer.write_text(json.dumps(manifest))
    before = source.read_bytes()
    parsed = ReportParser(build_report(tmp_path).read_text())
    assert len(parsed.images) == 1
    assert str(parsed.images[0]["src"]).endswith(".webp")
    with Image.open(tmp_path / str(parsed.images[0]["src"])) as small:
        assert small.width <= PREVIEW_SIZE[0] and small.height <= PREVIEW_SIZE[1]
    assert any(link.endswith("timecourse.png") for link in parsed.links)
    assert source.read_bytes() == before
    assert (
        tmp_path / str(parsed.images[0]["src"])
    ).stat().st_size < source.stat().st_size / 10


def test_preview_cache_reuses_images_after_relocation(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Cache hits need metadata only, including after moving a complete results tree."""
    root = tmp_path / "original"
    root.mkdir()
    source = root / "plot.png"
    with Image.new("RGB", (2000, 1000), "white") as image:
        image.save(source)
    first = preview_image(source, root)
    stat = first.path.stat()
    moved = tmp_path / "moved"
    shutil.copytree(root, moved)
    monkeypatch.setattr(
        "PIL.Image.open",
        lambda *args, **kwargs: pytest.fail("Cached previews must not decode images"),
    )
    again = preview_image(moved / "plot.png", moved)
    assert again.path.name == first.path.name
    assert (again.path.stat().st_size, again.path.stat().st_mtime_ns) == (
        stat.st_size,
        stat.st_mtime_ns,
    )


def test_changed_source_and_missing_preview_are_repaired(tmp_path: Path) -> None:
    """Invalidate changed sources and recover a missing preview without altering exports."""
    source = tmp_path / "plot.png"
    with Image.new("RGB", (2000, 1000), "white") as image:
        image.save(source)
    first = preview_image(source, tmp_path)
    first.path.unlink()
    assert preview_image(source, tmp_path).path.is_file()
    with Image.new("RGB", (2000, 1000), "red") as image:
        image.save(source)
    second = preview_image(source, tmp_path)
    assert second.path != first.path
    with Image.open(second.path) as small:
        pixel = small.getpixel((0, 0))
        assert isinstance(pixel, tuple)
        assert pixel[:3] == (255, 0, 0)


def test_preview_destinations_reject_symlinks(tmp_path: Path) -> None:
    """Never write cached previews through a link to unrelated output directories."""
    source = tmp_path / "plot.png"
    with Image.new("RGB", (24, 16), "white") as image:
        image.save(source)
    external = tmp_path / "unrelated"
    external.mkdir()
    (tmp_path / "report-previews").symlink_to(external, target_is_directory=True)
    with pytest.raises(ValueError, match="symbolic links"):
        preview_image(source, tmp_path)
    assert not list(external.iterdir())
