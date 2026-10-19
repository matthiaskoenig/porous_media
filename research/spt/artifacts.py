"""Publish reusable local SPT artifacts through validated atomic generations."""

import json
import os
import shutil
import tempfile
import uuid
from collections.abc import Callable, Generator, Mapping
from contextlib import contextmanager, suppress
from pathlib import Path


def validate_output(path: Path) -> None:
    """Reject destinations that traverse symbolic links before creating files."""
    if path.is_symlink() or path.parent.resolve() != path.parent.absolute():
        raise ValueError(f"Output must not follow symbolic links: {path}")


@contextmanager
def atomic_output(path: Path) -> Generator[Path]:
    """Publish a complete temporary file beside its destination or retain the old one."""
    validate_output(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        dir=path.parent, prefix=".partial-", suffix=path.suffix, delete=False
    ) as handle:
        temporary = Path(handle.name)
    try:
        yield temporary
        validate_output(path)
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def atomic_text(path: Path, text: str) -> None:
    """Replace a text file only after its UTF-8 contents are completely written."""
    with atomic_output(path) as temporary:
        temporary.write_text(text, encoding="utf-8")


def atomic_json(path: Path, value: object) -> None:
    """Replace a report only after a complete JSON file has been written."""
    atomic_text(path, json.dumps(value, indent=2, allow_nan=False) + "\n")


def fingerprint(path: Path) -> dict[str, int | str | None]:
    """Detect changes through metadata without hashing whole network files."""
    try:
        stat = path.stat()
    except FileNotFoundError:
        return {"path": str(path), "size": None, "mtime_ns": None}
    return {"path": str(path), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns}


def generation_path(directory: Path, manifest: Mapping[str, object]) -> Path:
    """Resolve a published UUID generation without following symbolic links."""
    name = manifest.get("generation")
    if not isinstance(name, str) or not name.startswith("generation-"):
        raise ValueError("Invalid artifact generation")
    uuid.UUID(name.removeprefix("generation-"))
    path = directory / name
    if path.is_symlink() or not path.is_dir():
        raise ValueError("Artifact generation missing or a symbolic link")
    return path


def records_current(records: object, root: Path | None = None) -> bool:
    """Check metadata records, optionally restricting them to a local output tree."""
    if not isinstance(records, list) or not records:
        return False
    for record in records:
        if not isinstance(record, dict) or not isinstance(record.get("path"), str):
            return False
        path = Path(record["path"])
        if root is not None and (
            not path.resolve().is_relative_to(root) or path.is_symlink()
        ):
            return False
        if fingerprint(path) != record:
            return False
    return True


@contextmanager
def run_lock(directory: Path) -> Generator[None]:
    """Hold an OS lock until artifact access completes, including after failure."""
    if (
        directory.resolve() != directory.absolute()
        or (directory / ".lock").is_symlink()
    ):
        raise ValueError(
            "Artifact directories and locks must not follow symbolic links"
        )
    directory.mkdir(parents=True, exist_ok=True)
    with (directory / ".lock").open("a+b") as handle:
        if os.name == "nt":
            import msvcrt

            handle.write(b"\0")
            handle.flush()
            handle.seek(0)
            msvcrt.locking(handle.fileno(), msvcrt.LK_NBLCK, 1)
        else:
            import fcntl

            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        try:
            yield
        finally:
            if os.name == "nt":
                handle.seek(0)
                msvcrt.locking(handle.fileno(), msvcrt.LK_UNLCK, 1)
            else:
                fcntl.flock(handle, fcntl.LOCK_UN)


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
        generation = generation_path(directory, manifest)
        if not records_current(manifest.get("sources"), root) or not records_current(
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
    metadata: dict[str, object] | None = None,
    *,
    catalog: Callable[[Path], Mapping[str, object]] | None = None,
) -> Path:
    """Publish a complete generation and retain the previous result on failure."""
    generation = directory / f"generation-{uuid.uuid4()}"
    old = None
    with suppress(OSError, ValueError):
        old = json.loads((directory / "manifest.json").read_text())
    staging.rename(generation)
    try:
        manifest: dict[str, object] = {
            "configuration": configuration,
            "generation": generation.name,
            "sources": records,
            "outputs": [
                fingerprint(path)
                for path in sorted(generation.rglob("*"))
                if path.is_file()
            ],
        }
        if metadata is not None:
            manifest["metadata"] = metadata
        if catalog is not None:
            manifest["catalog"] = dict(catalog(generation))
        if records != [fingerprint(path) for path in paths]:
            raise RuntimeError("Inputs changed before artifact publication")
        atomic_json(directory / "manifest.json", manifest)
    except BaseException:
        shutil.rmtree(generation)
        raise
    if isinstance(old, dict):
        with suppress(OSError, ValueError):
            previous = generation_path(directory, old)
            if previous != generation:
                shutil.rmtree(previous)
    return generation
