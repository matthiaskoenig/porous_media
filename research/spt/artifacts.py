"""Publish reusable local SPT artifacts through validated atomic generations."""

import json
import shutil
import uuid
from contextlib import suppress
from pathlib import Path

from research.spt.cache import _generation, _unchanged, fingerprint
from research.spt.cli import atomic_json


def current_manifest(
    directory: Path, configuration: dict[str, object], root: Path
) -> dict[str, object] | None:
    """Check completed artifacts without opening geometry or numerical arrays."""
    try:
        manifest = json.loads((directory / "manifest.json").read_text())
        if (
            not isinstance(manifest, dict)
            or manifest.get("configuration") != configuration
        ):
            return None
        generation = _generation(directory, manifest)
        if not _unchanged(manifest.get("sources"), root) or not _unchanged(
            manifest.get("outputs"), generation
        ):
            return None
        return manifest
    except (OSError, ValueError):
        return None


def publish(
    directory: Path,
    staging: Path,
    configuration: dict[str, object],
    paths: tuple[Path, ...],
    records: list[dict[str, int | str | None]],
    metadata: dict[str, object],
) -> Path:
    """Publish a complete generation and retain the previous result on failure."""
    generation = directory / f"generation-{uuid.uuid4()}"
    old = None
    with suppress(OSError, ValueError):
        old = json.loads((directory / "manifest.json").read_text())
    staging.rename(generation)
    try:
        manifest = {
            "configuration": configuration,
            "generation": generation.name,
            "sources": records,
            "outputs": [
                fingerprint(path)
                for path in sorted(generation.rglob("*"))
                if path.is_file()
            ],
            "metadata": metadata,
        }
        if records != [fingerprint(path) for path in paths]:
            raise RuntimeError("Inputs changed before artifact publication")
        atomic_json(directory / "manifest.json", manifest)
    except BaseException:
        shutil.rmtree(generation)
        raise
    if isinstance(old, dict):
        with suppress(OSError, ValueError):
            previous = _generation(directory, old)
            if previous != generation:
                shutil.rmtree(previous)
    return generation
