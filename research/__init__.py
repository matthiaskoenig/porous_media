"""Repository-only research workflows, separate from the distributed library."""

import os
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = Path(os.environ.get("POROUS_MEDIA_DATA_DIR", BASE_DIR / "data")).expanduser()
RESULTS_DIR = Path(
    os.environ.get("POROUS_MEDIA_RESULTS_DIR", BASE_DIR / "results")
).expanduser()
