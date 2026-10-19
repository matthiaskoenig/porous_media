"""Validate path components used for generated files."""

from ntpath import isreserved


def filename_component(value: str) -> str:
    """Require one portable filename component without directory traversal."""
    if (
        not value
        or value in {".", ".."}
        or any(character in value for character in "/\\:")
        or isreserved(value)
    ):
        raise ValueError(f"Expected a single portable filename component: {value!r}")
    return value
