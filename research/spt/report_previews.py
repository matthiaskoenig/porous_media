"""Cache small report previews independently of the original scientific exports."""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path

from research.spt.artifacts import atomic_json, atomic_output, validate_output

PREVIEW_SIZE = (768, 768)
PREVIEW_VERSION = 1


@dataclass(frozen=True)
class Preview:
    """Locate a cached preview and reserve its intrinsic dimensions in HTML."""

    path: Path
    width: int
    height: int


def preview_image(source: Path, root: Path) -> Preview:
    """Read a PNG only when its small, lossless WebP preview needs regeneration.

    Relative source paths and file metadata keep cache hits valid after relocation.
    Numerical arrays and raw source histories are never opened by this helper.
    """
    root = root.resolve()
    source = source.absolute()
    if not source.is_relative_to(root) or source.resolve() != source:
        raise ValueError("Report preview sources must stay inside the results root")
    stat = source.stat()
    configuration = {
        "version": PREVIEW_VERSION,
        "source": source.relative_to(root).as_posix(),
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
        "dimensions": PREVIEW_SIZE,
    }
    key = hashlib.sha256(
        json.dumps(configuration, sort_keys=True).encode()
    ).hexdigest()[:24]
    target = root / "report-previews" / f"{key}.webp"
    metadata = target.with_suffix(".json")
    validate_output(target)
    validate_output(metadata)
    if metadata.is_file() and target.is_file():
        try:
            cached = json.loads(metadata.read_text(encoding="utf-8"))
            image_stat = target.stat()
            if (cached["size"], cached["mtime_ns"]) == (
                image_stat.st_size,
                image_stat.st_mtime_ns,
            ):
                return Preview(target, int(cached["width"]), int(cached["height"]))
        except (OSError, ValueError, KeyError, TypeError):
            pass
    from PIL import Image

    with Image.open(source) as original:
        original.thumbnail(PREVIEW_SIZE, Image.Resampling.LANCZOS)
        with original.convert("RGBA" if "A" in original.getbands() else "RGB") as image:
            width, height = image.size
            with atomic_output(target) as temporary:
                image.save(temporary, format="WEBP", lossless=True, method=4)
    image_stat = target.stat()
    atomic_json(
        metadata,
        {
            "width": width,
            "height": height,
            "size": image_stat.st_size,
            "mtime_ns": image_stat.st_mtime_ns,
        },
    )
    return Preview(target, width, height)
