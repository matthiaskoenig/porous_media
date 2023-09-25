"""porous_media - Python utilities to work with porous media simulations."""

from pathlib import Path

__author__ = "Matthias König"
__version__ = "0.1.3"

program_name: str = "porous_media"

# data for analysis (not distributed)
BASE_DIR = Path(__file__).parent.parent.parent
DATA_DIR = BASE_DIR / "data"

# resources (distributed)
RESOURCES_DIR = Path(__file__).parent / "resources"

EXAMPLE_VTK = DATA_DIR / "example.vtk"
EXAMPLE_SIMULATION_XLSX = DATA_DIR / "simulation.xlsx"

# results directory
RESULTS_DIR = BASE_DIR / "results"
