"""Annotate existing manuscript curves and assemble the convergence supplement."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import numpy as np
from matplotlib.figure import Figure
from matplotlib.lines import Line2D

from research.spt.analysis_plots import LABELS, _scale
from research.spt.fields import LAYERS, SUMMARY_FIELDS, Layer
from research.spt.plot_data import PlotSource
from research.spt.plot_render import (
    RenderOptions,
    Samples,
    figure,
    label_sectors,
)
from research.spt.plot_style import convergence_color, simulation_color
from research.spt.publication_data import (
    CONVERGENCE_TIME,
    PATTERN_NAMES,
    PublicationData,
)
from research.spt.publication_geometry import (
    FIELD_TITLES,
    map_panel,
    raster_panel,
    save_publication,
)
from research.spt.publication_layout import resize_publication


def annotate_timecourses(
    art: Figure, axes: np.ndarray, *, width_mm: float | None = None
) -> None:
    """Add the manuscript's field headers, pattern names and dose-shade legend."""
    for column, title in enumerate(FIELD_TITLES):
        axes[0, column].set_title(title, fontsize=14, fontweight="bold", pad=12)
    for row, name in enumerate(PATTERN_NAMES):
        axis = axes[row, 0]
        axis.annotate(
            name.replace(" ", "\n", 1),
            xy=(0, 0.5),
            xycoords="axes fraction",
            xytext=(-65, 0),
            textcoords="offset points",
            rotation=90,
            ha="center",
            va="center",
            fontsize=12,
            fontweight="bold",
        )
        axis.annotate(
            chr(65 + row),
            xy=(0, 1),
            xycoords="axes fraction",
            xytext=(-76, 0),
            textcoords="offset points",
            ha="right",
            va="top",
            fontsize=16,
            fontweight="bold",
        )
    handles = [
        Line2D(
            [],
            [],
            color=simulation_color(0, dose),
            marker="o",
            markeredgecolor="black",
            label=str(dose),
        )
        for dose in range(1, 9)
    ]
    art.legend(
        handles=handles,
        title="Substrate inflow: dose index (low to high)",
        loc="outside upper center",
        ncol=8,
        frameon=False,
    )
    resize_publication(art, width_mm)


def render_convergence(
    data: PublicationData,
    sources: list[PlotSource],
    directory: Path,
    options: RenderOptions,
    *,
    comparison_time: float = CONVERGENCE_TIME,
    width_mm: float | None = None,
) -> None:
    """Compare meshes, common-time mixed fields and full native lobule curves."""
    mixed_samples = data.mixed_snapshot(comparison_time)
    mixed_sources = [
        source for source in sources if not source.run.run_id.endswith("50000")
    ]
    ranges = data.ranges(mixed_sources)
    art = figure((19, 15))
    outer = art.add_gridspec(4, 1, height_ratios=(2.4, 4.3, 3, 3), hspace=0.12)
    mesh_grid = outer[0].subgridspec(1, len(sources) + 1)
    common_end = min(source.end for source in sources)
    datasets = []
    try:
        for index, source in enumerate(sources):
            mesh = source.run.run_id.removeprefix("convergence_")
            lobulus = data.mesh(source.run.study, source.run.run_id)
            selected = lobulus.sectors[lobulus.cell_ids] == 0
            samples = Samples(
                lobulus.faces[selected],
                lobulus.cell_ids[selected],
                lobulus.values,
                lobulus.sectors,
                lobulus.times,
            )
            axis = art.add_subplot(mesh_grid[0, index])
            raster_panel(
                axis,
                samples,
                Layer("mesh", "Mesh", edge_alpha=0.45),
                np.zeros(len(samples.sectors)),
                (0, 1),
                side=2.3,
                dpi=options.dpi,
                mesh=True,
            )
            axis.set_title(
                f"{mesh}\n{np.count_nonzero(lobulus.sectors == 0):,} sixth cells",
                fontsize=10,
                fontweight="bold",
            )
            if index == 0:
                axis.text(
                    -0.2,
                    1.25,
                    "A",
                    transform=axis.transAxes,
                    fontsize=18,
                    fontweight="bold",
                )
            datasets.append(data.summary(source, native=True))
        mixed = data.mesh("convergence_combinations", "convergence_mixed")
        axis = art.add_subplot(mesh_grid[0, -1])
        raster_panel(
            axis,
            mixed,
            Layer("mesh", "Mesh", edge_alpha=0.45),
            np.zeros(len(mixed.sectors)),
            (0, 1),
            side=2.3,
            dpi=options.dpi,
            mesh=True,
        )
        label_sectors(axis, mixed)
        axis.set_title("Mixed sixths\n06250 to 37250", fontsize=10, fontweight="bold")
        maps = outer[1].subgridspec(1, len(LAYERS))
        for column, original in enumerate(LAYERS):
            layer = replace(
                original,
                show_edges=False,
                outline=True,
                limits=(0, 1) if original.sid == "rr_necrosis" else None,
            )
            axis = art.add_subplot(maps[0, column])
            map_panel(
                art,
                axis,
                mixed_samples,
                layer,
                mixed_samples.values[0, column],
                ranges,
                0.8,
                options.dpi,
                sectors=True,
            )
            axis.set_title(
                f"{FIELD_TITLES[column].replace(chr(10), ' ')} at {comparison_time / 3600:g} h",
                fontsize=12,
                fontweight="bold",
            )
            if column == 0:
                label_sectors(axis, mixed_samples)
                axis.text(
                    -0.2,
                    1.15,
                    "B",
                    transform=axis.transAxes,
                    fontsize=18,
                    fontweight="bold",
                )
        handles = []
        for row, symlog in enumerate((False, True), start=2):
            grid = outer[row].subgridspec(1, len(SUMMARY_FIELDS))
            for column, field in enumerate(SUMMARY_FIELDS):
                axis = art.add_subplot(grid[0, column])
                for index, (source, dataset) in enumerate(
                    zip(sources, datasets, strict=True)
                ):
                    mesh = source.run.run_id.removeprefix("convergence_")
                    color = convergence_color(mesh, index, len(sources))
                    time = dataset.time.values / 3600
                    values = dataset[f"{field}_mean"].values * (
                        100 if field == "rr_necrosis" else 1
                    )
                    (line,) = axis.plot(
                        time, values, color=color, linewidth=1.5, label=mesh
                    )
                    if row == 2 and column == 0:
                        handles.append(line)
                    if source.end < max(item.end for item in sources):
                        axis.plot(time[-1], values[-1], "x", color=color, markersize=5)
                _scale(axis, field, symlog)
                axis.axvline(
                    common_end / 3600, color="darkgray", linestyle=":", linewidth=0.8
                )
                axis.axvline(
                    comparison_time / 3600,
                    color="darkgray",
                    linestyle="--",
                    linewidth=0.8,
                )
                axis.set_xlim(0, max(source.end for source in sources) / 3600)
                axis.set_xticks([0, 6, 12, 18, 24])
                axis.set_xlabel("Time [hr]", fontsize=11, fontweight="bold")
                axis.set_ylabel(LABELS[field], fontsize=11, fontweight="bold")
                axis.set_title(
                    FIELD_TITLES[column].replace("\n", " ")
                    + (" (symlog)" if symlog and field != "rr_necrosis" else ""),
                    fontsize=12,
                    fontweight="bold",
                )
                if column == 0:
                    axis.text(
                        -0.23,
                        1.15,
                        "C" if row == 2 else "D",
                        transform=axis.transAxes,
                        fontsize=18,
                        fontweight="bold",
                    )
        art.legend(
            handles=handles,
            loc="outside upper center",
            ncol=len(handles),
            title="Convergence mesh identifier (coarse to fine)",
            frameon=False,
        )
        short = "; ".join(
            f"{source.run.run_id.removeprefix('convergence_')}: {source.end / 3600:.2f} h"
            for source in sources
            if source.end < max(item.end for item in sources)
        )
        art.supxlabel(
            f"A: meshes; B: mixed readouts at {comparison_time / 3600:g} h (dashed line); "
            "C: linear; D: symlog axes (necrosis stays linear). "
            f"Dotted line: common interval ends at {common_end / 3600:.2f} h.\n"
            f"Crosses mark truncated native histories ({short}); no extrapolation. "
            "Mixed sectors 1-6: 06250, 09375, 12500, 18750, 25000, 37250.",
            fontsize=10,
        )
        save_publication(
            art, directory / "figure-s1-convergence", options, width_mm=width_mm
        )
    finally:
        for dataset in datasets:
            dataset.close()
        art.clear()
