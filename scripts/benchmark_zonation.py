"""Measure zonation on a reproducible synthetic mesh without rendering."""

import argparse
import json
import statistics
import time

import meshio
import numpy as np

from porous_media.mesh.mesh_zonation import ZonatedMesh


def benchmark(cells: int, repeat: int) -> dict[str, int | float]:
    """Return median wall time for a mesh with two boundary regions."""
    if cells < 4 or repeat < 1:
        raise ValueError("cells must be at least 4 and repeat must be positive")
    centers = np.linspace(0.0, 1.0, cells)
    points = np.zeros((2 * cells, 3))
    points[::2, 0] = centers - 1e-6
    points[1::2, 0] = centers + 1e-6
    labels = np.zeros(cells)
    labels[: cells // 4] = 1
    labels[-(cells // 4) :] = 2
    mesh = meshio.Mesh(
        points,
        [("line", np.arange(2 * cells).reshape(-1, 2))],
        cell_data={"cell_type": [labels], "element_volume_point_TPM": [np.ones(cells)]},
    )
    durations = []
    for _ in range(repeat):
        start = time.perf_counter()
        ZonatedMesh.create_zonated_mesh(mesh)
        durations.append(time.perf_counter() - start)
    return {
        "cells": cells,
        "repeat": repeat,
        "median_seconds": statistics.median(durations),
    }


def main() -> None:
    """Print benchmark results as JSON for local comparisons."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cells", type=int, default=2000)
    parser.add_argument("--repeat", type=int, default=3)
    args = parser.parse_args()
    print(json.dumps(benchmark(args.cells, args.repeat)))


if __name__ == "__main__":
    main()
