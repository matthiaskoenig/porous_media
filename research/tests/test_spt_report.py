"""Check complete, offline, relocatable SPT inspection reports with synthetic data."""

import csv
import json
import shutil
import uuid
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote, urlsplit

import pytest

from research.spt.report import build_report, main


class ReportParser(HTMLParser):
    """Collect cards, resources, and links through the browser's HTML syntax."""

    def __init__(self, document: str) -> None:
        """Parse a complete generated document."""
        super().__init__(convert_charrefs=True)
        self.cards: list[dict[str, str | None]] = []
        self.links: list[str] = []
        self.images: list[dict[str, str | None]] = []
        self.identifiers: list[str] = []
        self.feed(document)

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        """Track published figure cards and their relative resources."""
        values = dict(attrs)
        if values.get("id"):
            self.identifiers.append(str(values["id"]))
        if tag == "article" and values.get("class") == "figure-card":
            self.cards.append(values)
        if tag == "a" and values.get("href"):
            self.links.append(str(values["href"]))
        if tag == "img" and values.get("src"):
            self.images.append(values)


def _publish(
    root: Path,
    directory: str,
    names: tuple[str, ...],
    configuration: dict[str, object] | None = None,
) -> Path:
    path = root / directory
    generation = path / f"generation-{uuid.uuid4()}"
    generation.mkdir(parents=True)
    records = []
    for name in names:
        output = generation / name
        output.write_text("{}" if output.suffix == ".json" else "synthetic output")
        records.append(
            {
                "path": str(output),
                "size": output.stat().st_size,
                "mtime_ns": output.stat().st_mtime_ns,
            }
        )
    manifest = path / "manifest.json"
    manifest.write_text(
        json.dumps(
            {
                "generation": generation.name,
                "outputs": records,
                "configuration": configuration or {},
                # Deliberately unavailable. Reports must not follow source paths.
                "sources": [{"path": "/unavailable/network/source.h5"}],
            }
        )
    )
    return generation


@pytest.fixture
def published(tmp_path: Path) -> Path:
    """Create each report section without rendering or valid numerical files."""
    root = tmp_path / "results"
    _publish(
        root,
        "analysis/simulations_lobulus",
        ("timecourses.png", "timecourses.svg", "analysis.nc", "native.tsv"),
        {"study": "simulations_lobulus", "runs": ["sim001"]},
    )
    for geometry in ("sixth", "lobulus"):
        _publish(
            root,
            f"figures/simulations_{geometry}/sim001",
            ("timecourse.png", "timecourse.pdf", "provenance.json"),
            {"title": "sim001 | constant | dose 1"},
        )
    _publish(root, "overviews/simulations_lobulus", ("meshes.svg",))
    _publish(root, "simulations_combinations/xdmf/dose_01", ("dose_01.png",))
    _publish(
        root,
        "convergence_combinations/xdmf/convergence_mixed",
        ("convergence_mixed.png",),
    )
    _publish(
        root,
        "convergence-symlog/figures/convergence_sixth/convergence_06250",
        ("timecourse.png",),
        {"options": {"layers": [{"title": "Substrate", "scale": "symlog"}]}},
    )
    (root / "figures/boundary_color.png").write_text("condition figure")
    with (root / "catalog.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=(
                "run_id",
                "pattern_name",
                "boundary_flow_key",
                "complete",
                "tend_s",
            ),
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerow(
            {
                "run_id": "sim001",
                "pattern_name": "constant",
                "boundary_flow_key": "0.0",
                "complete": "True",
                "tend_s": "4",
            }
        )
        writer.writerow(
            {"run_id": "convergence_06250", "complete": "False", "tend_s": "2"}
        )
    return root


def test_report_covers_all_figures_and_exports(published: Path) -> None:
    """Group exports, preserve every section, and include data and quality details."""
    before = {
        path: (path.stat().st_size, path.stat().st_mtime_ns)
        for path in published.rglob("*")
        if path.is_file()
    }
    document = build_report(published).read_text()
    parsed = ReportParser(document)
    assert len(parsed.cards) == 8
    assert len({card["id"] for card in parsed.cards}) == 8
    assert len(set(parsed.identifiers)) == len(parsed.identifiers)
    assert {card["data-section"] for card in parsed.cards} == {
        "analysis",
        "geometry",
        "combinations",
        "overviews",
        "supplements",
        "conditions",
    }
    mixed = next(
        card
        for card in parsed.cards
        if card["data-title"] == "Six convergence resolutions"
    )
    assert mixed["data-geometry"] == "Combined"
    assert mixed["data-family"] == "Convergence"
    simulation = next(
        card for card in parsed.cards if card["data-section"] == "analysis"
    )
    assert json.loads(str(simulation["data-pattern"])) == ["constant"]
    assert json.loads(str(simulation["data-dose"])) == ["1"]
    assert all(image.get("loading") == "lazy" for image in parsed.images)
    assert not any(str(image["src"]).endswith(".svg") for image in parsed.images)
    assert any(link.endswith("analysis.nc") for link in parsed.links)
    assert any(link.endswith("native.tsv") for link in parsed.links)
    assert "convergence_06250 ends at 2 s" in document
    assert "omitted generations" not in document
    assert (published / "index.html").read_text() == document
    assert {
        path: (path.stat().st_size, path.stat().st_mtime_ns) for path in before
    } == before
    for link in parsed.links + [str(image["src"]) for image in parsed.images]:
        url = urlsplit(link)
        assert not url.scheme and not url.netloc
        if url.path:
            assert (published / unquote(url.path)).is_file()


def test_offline_relocation_and_stable_figure_links(
    published: Path, tmp_path: Path
) -> None:
    """Follow only local generation suffixes after moving a complete result tree."""
    original = ReportParser(build_report(published).read_text())
    moved = tmp_path / "moved results"
    shutil.copytree(published, moved)
    shutil.rmtree(published)
    document = build_report(moved).read_text()
    parsed = ReportParser(document)
    assert [card["id"] for card in parsed.cards] == [
        card["id"] for card in original.cards
    ]
    assert "omitted generations" not in document
    assert str(published) not in document
    assert "/unavailable/network" not in document
    _publish(moved, "analysis/convergence_lobulus", ("additional.png",))
    refreshed = ReportParser(build_report(moved).read_text())
    assert {card["id"] for card in parsed.cards} < {
        card["id"] for card in refreshed.cards
    }


@pytest.mark.parametrize("damage", ["missing", "changed", "traversal", "symlink"])
def test_invalid_generations_are_omitted(
    published: Path, tmp_path: Path, damage: str
) -> None:
    """Report omissions and retain valid sections without unsafe or broken links."""
    directory = published / "analysis/simulations_lobulus"
    manifest = directory / "manifest.json"
    value = json.loads(manifest.read_text())
    image = directory / value["generation"] / "timecourses.png"
    if damage == "missing":
        image.unlink()
    elif damage == "changed":
        image.write_text("modified")
    elif damage == "traversal":
        value["outputs"][0]["path"] = str(image.parent / ".." / "outside.png")
        manifest.write_text(json.dumps(value))
    else:
        external = tmp_path / "external.png"
        image.rename(external)
        image.symlink_to(external)
    document = build_report(published).read_text()
    parsed = ReportParser(document)
    assert len(parsed.cards) == 7
    assert "1 omitted generation requires attention" in document
    assert not any(link.endswith("timecourses.png") for link in parsed.links)


def test_escape_metadata_and_encode_filenames(tmp_path: Path) -> None:
    """Never turn figure metadata into markup, scripts, or remote resource URLs."""
    root = tmp_path / "report"
    title = '</script><script>alert("title")</script><img onerror="alert(1)">'
    _publish(
        root,
        "figures/simulations_sixth/sim001",
        ('space # & " figure.png',),
        {"title": title},
    )
    document = build_report(root).read_text()
    parsed = ReportParser(document)
    assert len(parsed.cards) == 1
    assert parsed.cards[0]["data-title"] == title
    assert title not in document
    assert len(parsed.images) == 1
    assert "%23%20%26%20%22" in str(parsed.images[0]["src"])
    assert (root / unquote(str(parsed.images[0]["src"]))).is_file()


def test_report_destinations_and_missing_results(tmp_path: Path) -> None:
    """Keep prior entry points on failure and never overwrite symlink targets."""
    missing = tmp_path / "missing"
    assert main(["--output-root", str(missing)]) == 1
    assert not missing.exists()
    root = tmp_path / "results"
    root.mkdir()
    old = root / "report.html"
    old.write_text("previous report")
    target = tmp_path / "unrelated.html"
    target.write_text("unrelated")
    (root / "index.html").symlink_to(target)
    with pytest.raises(ValueError, match="symbolic link"):
        build_report(root)
    assert old.read_text() == "previous report"
    assert target.read_text() == "unrelated"
    (root / "index.html").unlink()
    assert main(["--output-root", str(root)]) == 0
    assert "No published figures are available" in old.read_text()
