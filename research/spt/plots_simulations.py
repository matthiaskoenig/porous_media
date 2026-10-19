"""Analyze SPT simulations lobules with the shared reproducible workflow.

Run with ``--help`` for source selection, spatial profiles and serialized results.
"""

from research.spt.analysis import main

if __name__ == "__main__":
    raise SystemExit(main(study="simulations"))
