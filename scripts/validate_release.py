"""Validate a release version and its notes before building distributions."""

import argparse
from pathlib import Path

from porous_media import __version__


def validate_release(tag: str, notes_dir: Path) -> Path:
    """Require a matching version and nonempty UTF-8 release notes.

    Args:
        tag: Release tag without a leading ``v``.
        notes_dir: Directory containing versioned Markdown release notes.

    Returns:
        Path to the validated release notes.

    Raises:
        ValueError: The tag differs from the version or the notes are empty.
        OSError: The release notes cannot be read.
    """
    if tag != __version__:
        raise ValueError(f"Tag {tag!r} differs from package version {__version__!r}")
    notes = notes_dir / f"{tag}.md"
    if not notes.read_text(encoding="utf-8").strip():
        raise ValueError(f"Empty release notes: {notes}")
    return notes


def main() -> None:
    """Validate a prospective release without creating tags or publishing."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("tag", help="Release version without a leading v")
    parser.add_argument(
        "--notes-dir",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "docs" / "release-notes",
    )
    args = parser.parse_args()
    try:
        notes = validate_release(args.tag, args.notes_dir)
    except (OSError, ValueError) as error:
        parser.error(str(error))
    print(f"Validated release {args.tag}: {notes}")


if __name__ == "__main__":
    main()
