"""Exercise the release validation command without publishing artifacts."""

import subprocess
import sys
from pathlib import Path

import pytest

from porous_media import __version__

VALIDATOR = Path(__file__).resolve().parents[1] / "scripts" / "validate_release.py"


def run_validator(tag: str, notes_dir: Path) -> subprocess.CompletedProcess[str]:
    """Run the same validation command used by the build workflow."""
    return subprocess.run(
        [
            sys.executable,
            str(VALIDATOR),
            tag,
            "--notes-dir",
            str(notes_dir),
        ],
        cwd=notes_dir,
        capture_output=True,
        text=True,
        check=False,
    )


def test_release_validation_accepts_matching_version_and_notes(tmp_path: Path) -> None:
    """Validate Unicode notes from outside the repository working directory."""
    (tmp_path / f"{__version__}.md").write_text(
        f"# {__version__}\n\nChanges by König.\n", encoding="utf-8"
    )
    result = run_validator(__version__, tmp_path)
    assert result.returncode == 0, result.stderr
    assert f"Validated release {__version__}" in result.stdout


@pytest.mark.parametrize("tag", ["0.0.0", f"v{__version__}", "../outside"])
def test_release_validation_rejects_incorrect_tag(tmp_path: Path, tag: str) -> None:
    """Reject a mismatched version before accessing release notes."""
    result = run_validator(tag, tmp_path)
    assert result.returncode != 0
    assert f"Tag {tag!r} differs from package version {__version__!r}" in result.stderr


def test_release_validation_rejects_missing_notes(tmp_path: Path) -> None:
    """Fail the build when the versioned release notes are absent."""
    result = run_validator(__version__, tmp_path)
    assert result.returncode != 0
    assert f"{__version__}.md" in result.stderr


@pytest.mark.parametrize("content", ["", " \n\t"])
def test_release_validation_rejects_empty_notes(tmp_path: Path, content: str) -> None:
    """Require release notes to describe changes before publishing."""
    (tmp_path / f"{__version__}.md").write_text(content, encoding="utf-8")
    result = run_validator(__version__, tmp_path)
    assert result.returncode != 0
    assert "Empty release notes" in result.stderr
