"""Create SPT simulation figures for sixth, lobulus and combined geometries.

Use ``--help`` for field/time selection, vector formats and opt-in videos.
"""

from research.spt.plot_geometry import main
from research.spt.plot_legacy import necrosis_plots, visualize_spt_2d

__all__ = ["necrosis_plots", "visualize_spt_2d"]

if __name__ == "__main__":
    raise SystemExit(main(study="simulations"))
