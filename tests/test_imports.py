"""Keep library imports independent of research configuration and application state."""

import subprocess
import sys
from pathlib import Path

import pytest


@pytest.mark.parametrize(
    "code",
    [
        "import sys; import numpy as np; np.random.seed(123); "
        "expected = np.random.RandomState(123).random_sample(); "
        "import porous_media.mesh.mesh_zonation; "
        "assert np.random.random_sample() == expected; "
        "assert 'pyvista' not in sys.modules",
        "import pyvista as pv; pv.global_theme.background = 'red'; "
        "pv.global_theme.font.size = 17; "
        "import porous_media.visualization.pyvista_visualization; "
        "assert pv.global_theme.background == pv.Color('red'); "
        "assert pv.global_theme.font.size == 17",
        "import importlib.util; import porous_media; "
        "assert importlib.util.find_spec('porous_media.analyses') is None; "
        "assert not hasattr(porous_media, 'DATA_DIR'); "
        "assert not hasattr(porous_media, 'RESULTS_DIR')",
        "import sys; hook = sys.displayhook; import porous_media.console; "
        "assert sys.displayhook is hook; "
        "assert not porous_media.console.console.record",
        "import matplotlib; matplotlib.rcParams['font.size'] = 42; "
        "from porous_media.visualization.plots import plot_style; "
        "assert matplotlib.rcParams['font.size'] == 42; "
        "context = plot_style(); context.__enter__(); "
        "assert matplotlib.rcParams['font.size'] == 13; "
        "context.__exit__(None, None, None); "
        "assert matplotlib.rcParams['font.size'] == 42",
    ],
)
def test_library_import_is_quiet_and_preserves_state(tmp_path: Path, code: str) -> None:
    """Import from an arbitrary working directory without research data or side effects."""
    result = subprocess.run(
        [sys.executable, "-c", code],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout == ""
