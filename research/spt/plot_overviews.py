"""Compose final-state and distinct-mesh comparisons directly from sampled arrays."""

from __future__ import annotations

import json
import math
import tempfile
from dataclasses import asdict
from pathlib import Path

import matplotlib as mpl
import numpy as np

from research.spt.artifacts import (
    current_manifest,
    fingerprint,
    generation_path,
    publish,
    records_current,
    run_lock,
)
from research.spt.conditions import pattern_idx2name
from research.spt.fields import Layer
from research.spt.plot_data import PlotSource, load_heights
from research.spt.plot_mesh import OUTLINE_VERSION
from research.spt.plot_render import (
    RENDER_VERSION,
    RenderOptions,
    Samples,
    display_limits,
    figure,
    label_sectors,
    panel,
    save_figure,
    scalar_bar,
)
from research.spt.plot_view import GEOMETRY_VIEW

OVERVIEW_VERSION = 4


def _label(source: PlotSource) -> str:
    if source.family == "convergence":
        return source.run.run_id.removeprefix("convergence_")
    return f"Dose {source.dose + 1}" if source.dose is not None else source.run.run_id


def _layout(
    sources: list[PlotSource],
) -> tuple[int, int, list[tuple[int, int]], list[str], list[str]]:
    if (
        sources[0].family == "simulations"
        and sources[0].geometry != "combined"
        and all(source.run.condition and source.dose is not None for source in sources)
    ):
        patterns = sorted(
            {
                source.run.condition.pattern_key
                for source in sources
                if source.run.condition
            }
        )
        doses = sorted({source.dose for source in sources if source.dose is not None})
        positions = [
            (patterns.index(source.run.condition.pattern_key), doses.index(source.dose))
            for source in sources
            if source.run.condition and source.dose is not None
        ]
        return (
            len(patterns),
            len(doses),
            positions,
            [pattern_idx2name[key] for key in patterns],
            [f"Dose {key + 1}" for key in doses],
        )
    columns = min(4, len(sources))
    return (
        math.ceil(len(sources) / columns),
        columns,
        [divmod(index, columns) for index in range(len(sources))],
        [],
        [],
    )


def render_overview(
    sources: list[PlotSource],
    samples: list[Samples],
    directory: Path,
    options: RenderOptions,
    ranges: dict[str, tuple[float, float]],
) -> None:
    """Render final states in condition order and each distinct geometry only once."""
    preferred = "rr_necrosis" if sources[0].family == "simulations" else "rr_(T)"
    layer = next(
        (layer for layer in options.layers if layer.sid == preferred),
        options.layers[-1],
    )
    field = options.layers.index(layer)
    rows, columns, positions, row_labels, column_labels = _layout(sources)
    art = figure((4.0 * columns, 4.0 * rows))
    axes = art.subplots(rows, columns, squeeze=False)
    try:
        for axis in axes.flat:
            axis.set_axis_off()
        for source, data, (row, column) in zip(
            sources, samples, positions, strict=True
        ):
            axis = axes[row, column]
            collection = panel(
                axis,
                data,
                layer,
                data.frame(options.snapshots[-1])[field],
                display_limits(layer, ranges),
                sectors=source.geometry == "combined",
            )
            scalar_bar(art, axis, collection, layer)
            if not column_labels:
                axis.set_title(_label(source), fontsize=11)
            elif row == 0:
                axis.set_title(column_labels[column], fontsize=11)
            if row_labels and column == 0:
                axis.text(
                    -0.03,
                    0.5,
                    row_labels[row],
                    transform=axis.transAxes,
                    ha="right",
                    va="center",
                    rotation=90,
                    fontsize=10,
                )
        if sources[0].sector_labels:
            label_sectors(axes.flat[0], samples[0])
            art.supxlabel(
                " | ".join(
                    f"{index + 1}: {name}"
                    for index, name in enumerate(sources[0].sector_labels)
                ),
                fontsize=9,
            )
        art.suptitle(
            f"{sources[0].family.title()} | {sources[0].geometry} | {layer.title} at {options.snapshots[-1] / 3600:g} h",
            fontsize=13,
        )
        save_figure(art, directory / "final-state", options)
    finally:
        art.clear()

    unique = []
    seen = set()
    for index, data in enumerate(samples):
        key = (data.faces.shape, data.faces.tobytes(), data.cell_ids.tobytes())
        if key not in seen:
            seen.add(key)
            unique.append(index)
    columns = min(4, len(unique))
    art = figure((3.5 * columns, 3.5 * math.ceil(len(unique) / columns) + 0.6))
    axes = art.subplots(math.ceil(len(unique) / columns), columns, squeeze=False)
    try:
        for axis in axes.flat:
            axis.set_axis_off()
        for axis, index in zip(axes.flat, unique, strict=False):
            data, source = samples[index], sources[index]
            values = np.maximum(data.sectors, 0)
            panel(
                axis,
                data,
                Layer("mesh", "Mesh", edge_alpha=0.45),
                values,
                (-0.5, 5.5),
                mesh=True,
            )
            if source.sector_labels:
                label_sectors(axis, data)
            axis.set_title(
                f"{source.run.run_id}\n{data.values.shape[2]:,} cells", fontsize=11
            )
        art.suptitle(
            f"{sources[0].family.title()} | {sources[0].geometry} | {len(unique)} distinct mesh(es)",
            fontsize=13,
        )
        if sources[0].sector_labels:
            art.supxlabel(
                " | ".join(
                    f"{index + 1}: {name}"
                    for index, name in enumerate(sources[0].sector_labels)
                ),
                fontsize=9,
            )
        save_figure(art, directory / "meshes", options)
    finally:
        art.clear()


def overview_figures(
    sources: list[PlotSource],
    prepared: list[dict[str, object]],
    output_root: Path,
    options: RenderOptions,
    ranges: dict[str, tuple[float, float]],
    *,
    force: bool = False,
) -> list[dict[str, str]]:
    """Publish one cached overview set per study and geometry."""
    source = sources[0]
    directory = source.directory(output_root, "overviews").parent
    paths = tuple(
        generation_path(item.directory(output_root, "plot-data"), manifest)
        / "samples.npz"
        for item, manifest in zip(sources, prepared, strict=True)
    )
    configuration: dict[str, object] = {
        "version": RENDER_VERSION,
        "geometry_view": GEOMETRY_VIEW,
        "overview_version": OVERVIEW_VERSION,
        "sources": [item.title for item in sources],
        "labels": list(source.sector_labels),
        "options": json.loads(json.dumps(asdict(options))),
        "ranges": {name: list(value) for name, value in ranges.items()},
    }
    if any(layer.outline for layer in options.layers):
        configuration["outline_version"] = OUTLINE_VERSION
    with run_lock(directory):
        if any(
            not records_current(manifest.get("sources"), item.run.root)
            for item, manifest in zip(sources, prepared, strict=True)
        ):
            raise ValueError("Overview samples are stale; prepare them again")
        cached = current_manifest(directory, configuration, output_root)
        if cached and not force:
            generation = generation_path(directory, cached)
        else:
            before = [fingerprint(path) for path in paths]
            samples = [
                Samples.load(
                    path,
                    heights=load_heights(
                        item.directory(output_root, "plot-depth"),
                        manifest,
                        item.run.root,
                        sample_path=path,
                    ),
                )
                for path, item, manifest in zip(paths, sources, prepared, strict=True)
            ]
            with tempfile.TemporaryDirectory(
                prefix=".partial-", dir=directory
            ) as temporary:
                staging = Path(temporary)
                with mpl.rc_context({"svg.fonttype": "none", "pdf.fonttype": 42}):
                    render_overview(sources, samples, staging, options, ranges)
                if any(
                    not records_current(manifest.get("sources"), item.run.root)
                    for item, manifest in zip(sources, prepared, strict=True)
                ):
                    raise RuntimeError("Overview inputs changed during rendering")
                generation = publish(
                    directory, staging, configuration, paths, before, {}
                )
    result = []
    for name in ("final-state", "meshes"):
        row = {"study": source.run.study, "run_id": name, "geometry": source.geometry}
        row.update(
            {
                extension: str(generation / f"{name}.{extension}")
                for extension in options.formats
            }
        )
        result.append(row)
    return result
