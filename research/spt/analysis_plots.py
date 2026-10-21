"""Compose SPT time, position, dose and mesh comparisons from serialized results."""

from __future__ import annotations

import html
import tempfile
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from typing import cast

import matplotlib as mpl
import numpy as np
import pandas as pd
import xarray as xr
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.ticker import FixedLocator, SymmetricalLogLocator

from research.spt.analysis_data import load_analysis, weighted_statistics
from research.spt.artifacts import (
    atomic_json,
    atomic_text,
    current_manifest,
    fingerprint,
    generation_path,
    publish,
    run_lock,
)
from research.spt.cache import open_processed_cells
from research.spt.conditions import pattern_idx2name
from research.spt.fields import DIAGNOSTICS, POSITION, VOLUME
from research.spt.plot_data import PlotSource
from research.spt.plot_render import RenderOptions, figure, save_figure
from research.spt.plot_style import (
    CONVERGENCE_ORDER,
    LABEL_STYLE,
    convergence_color,
    simulation_color,
)

ANALYSIS_PLOT_VERSION = 5
LABELS = {
    "rr_(S_ext)": "Substrate [mM]",
    "rr_(P_ext)": "Product [mM]",
    "rr_protein": "Protein [-]",
    "rr_(T)": "Toxic compound [mM]",
    "rr_necrosis": "Necrosis [%]",
    VOLUME: "Mesh cell volume [nl]",
    "rr_Vext": "Fluid volume [nl]",
    "rr_Vli": "Solid volume [nl]",
    "pressure": "Pressure [?]",
}
# These are the unit conversions declared in the manuscript's diagnostic script.
DIAGNOSTIC_FACTORS = {VOLUME: 1e12, "rr_Vext": 1e9, "rr_Vli": 1e9}


@dataclass(frozen=True)
class PlotOptions:
    """Select manuscript figures and optional population SD or symlog comparisons."""

    formats: tuple[str, ...] = ("svg", "pdf", "png")
    dpi: int = 300
    spread: bool | None = None
    convergence_symlog: bool = True

    def validate(self) -> None:
        """Reuse the shared format validation before writing any output."""
        RenderOptions(formats=self.formats, dpi=self.dpi).validate()


def comparisons(
    aligned: xr.Dataset, profiles: xr.Dataset, sources: list[PlotSource]
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compare meshes to the finest available mesh at identical times and positions."""
    reference = int(np.argmax(aligned.num_cells.values))
    ref_id = sources[reference].run.run_id
    delta = abs(aligned["mean"] - aligned["mean"].isel(run=reference))
    denominator = abs(aligned["mean"].isel(run=reference))
    differences = xr.Dataset(
        {
            "absolute_difference": delta,
            "relative_difference": delta / denominator.where(denominator > 0),
        }
    )
    means = differences.to_dataframe().reset_index()
    means["reference_run"] = ref_id
    rows = []
    for source in sources:
        run = source.run.run_id
        for time in profiles.time.values:
            for field in aligned.field.values:
                ref = profiles.sel(run=ref_id, time=time, field=field)
                current = profiles.sel(run=run, time=time, field=field)
                valid = (
                    np.isfinite(ref["mean"])
                    & np.isfinite(current["mean"])
                    & (ref.volume > 0)
                )
                matched = float(ref.volume.where(valid, 0).sum())
                total = float(ref.volume.sum())
                rmse = (
                    float(
                        np.sqrt(
                            (((current["mean"] - ref["mean"]) ** 2) * ref.volume)
                            .where(valid, 0)
                            .sum()
                            / matched
                        )
                    )
                    if matched > 0
                    else float("nan")
                )
                rows.append(
                    {
                        "run": run,
                        "time": float(time),
                        "field": str(field),
                        "reference_run": ref_id,
                        "profile_rmse": rmse,
                        "reference_volume_coverage": matched / total,
                    }
                )
    return means, pd.DataFrame(rows)


def _factor(field: str) -> float:
    return 100.0 if field == "rr_necrosis" else 1.0


def _style(source: PlotSource, index: int, count: int) -> tuple[str, str]:
    if source.run.condition:
        condition = source.run.condition
        return source.run.run_id, simulation_color(
            condition.pattern_key, condition.boundary_flow_key + 1
        )
    mesh = source.run.run_id.removeprefix("convergence_")
    return mesh, convergence_color(mesh, index, count)


def _scale(axis: Axes, field: str, symlog: bool) -> None:
    values = (
        np.concatenate(
            [
                np.asarray(line.get_ydata(), dtype=float).reshape(-1)
                for line in axis.lines
            ]
        )
        if axis.lines
        else np.array([0.0])
    )
    finite = values[np.isfinite(values)]
    if len(finite) and finite.min() >= 0:
        axis.set_ylim(bottom=0)
    if field == "rr_necrosis":
        axis.set_ylim(-2, 102)
    elif symlog:
        extent = max(
            (
                float(np.nanmax(abs(np.asarray(line.get_ydata(), dtype=float))))
                for line in axis.lines
                if np.asarray(line.get_ydata()).size
            ),
            default=1.0,
        )
        threshold = extent * 1e-6 if extent > 0 else 1e-6
        axis.set_yscale("symlog", linthresh=threshold)
        locator = SymmetricalLogLocator(base=10, linthresh=threshold)
        locator.set_params(numticks=5)
        low, high = axis.get_ylim()
        axis.yaxis.set_major_locator(
            FixedLocator(
                [
                    value
                    for value in locator.tick_values(low, high)
                    if low <= value <= high and (value == 0 or abs(value) >= threshold)
                ]
            )
        )
    else:
        axis.ticklabel_format(
            axis="y", style="sci", scilimits=(-2, 3), useMathText=True
        )


def _curves(
    sources: list[PlotSource],
    datasets: list[xr.Dataset],
    directory: Path,
    name: str,
    options: PlotOptions,
    fields: tuple[str, ...],
    *,
    spatial_time: float | None = None,
    symlog: bool = False,
    decorate: Callable[[Figure, np.ndarray], None] | None = None,
) -> None:
    """Keep manuscript panels, matched-time error bars and individual cell scatters."""
    simulation = sources[0].family == "simulations"
    patterns = (
        sorted(
            {
                source.run.condition.pattern_key
                for source in sources
                if source.run.condition
            }
        )
        if simulation
        else [0]
    )
    spread = simulation if options.spread is None else options.spread
    art = figure((2.7 * len(fields), 2.5 * len(patterns)))
    axes = art.subplots(len(patterns), len(fields), squeeze=False)
    maxima = dict.fromkeys(fields, 0.0)
    legend = {}
    try:
        for index, (source, dataset) in enumerate(zip(sources, datasets, strict=True)):
            row = (
                patterns.index(source.run.condition.pattern_key)
                if source.run.condition
                else 0
            )
            label, color = _style(source, index, len(sources))
            for column, field in enumerate(fields):
                axis = axes[row, column]
                factor = _factor(field)
                if spatial_time is None:
                    selected = dataset.sel(field=field)
                    x = selected.time.values / 3600
                    y, std = (
                        selected["mean"].values * factor,
                        selected["std"].values * factor,
                    )
                    artist = axis.errorbar(
                        x,
                        y,
                        yerr=std if spread else None,
                        color=color,
                        label=label,
                        linestyle="-",
                        linewidth=1.5,
                        marker="o" if simulation else "None",
                        markersize=6 if simulation else 3,
                        markeredgecolor="black",
                        markeredgewidth=1 if simulation else 0.5,
                    )
                    upper = y + std
                else:
                    selected = dataset.sel(time=spatial_time)
                    x, y = selected[POSITION].values, selected[field].values * factor
                    (artist,) = axis.plot(
                        x,
                        y,
                        color=color,
                        label=label,
                        linestyle="",
                        marker="o",
                        markersize=6 if simulation else 3,
                        markeredgecolor="black",
                        markeredgewidth=1 if simulation else 0.5,
                    )
                    upper = y
                finite = upper[np.isfinite(upper)]
                if len(finite):
                    maxima[field] = max(maxima[field], float(finite.max()))
                legend[label] = artist
        for row in range(len(patterns)):
            for column, field in enumerate(fields):
                axis = axes[row, column]
                axis.set_ylabel(LABELS[field], **LABEL_STYLE)
                maximum = maxima[field] or 1.0
                if not simulation and spatial_time is None:
                    axis.set_ylim(bottom=0)
                else:
                    axis.set_ylim(-0.05 * maximum if simulation else 0, 1.05 * maximum)
                if field == "rr_necrosis" and spatial_time is None and simulation:
                    axis.set_ylim(0, 105)
                if symlog:
                    _scale(axis, field, True)
                bottom = row == len(patterns) - 1
                if spatial_time is None:
                    if simulation:
                        axis.set_xticks([0, 10, 20])
                    else:
                        axis.axvline(
                            1.5, color="darkgray", linestyle="--", linewidth=0.5
                        )
                    if bottom:
                        axis.set_xlabel("Time [hr]", **LABEL_STYLE)
                else:
                    axis.set_xlim(0, 1)
                    axis.set_xticks([0, 0.5, 1], labels=["PP", "", "PV"])
                    if bottom:
                        axis.set_xlabel("Position [-]", **LABEL_STYLE)
                axis.tick_params(axis="x", labelbottom=bottom)
        if not simulation and spatial_time is None:
            labels = sorted(
                legend,
                key=lambda label: (
                    CONVERGENCE_ORDER.index(label)
                    if label in CONVERGENCE_ORDER
                    else len(CONVERGENCE_ORDER)
                ),
            )
            art.legend(
                [legend[label] for label in labels],
                labels,
                loc="lower center",
                bbox_to_anchor=(0.5, 1.0),
                ncol=len(legend),
                frameon=True,
            )
        if decorate is not None:
            decorate(art, axes)
        save_figure(
            art,
            directory / name,
            RenderOptions(formats=options.formats, dpi=options.dpi),
        )
    finally:
        art.clear()


def _initial_cells(
    source: PlotSource, cells: xr.Dataset, data_root: Path
) -> xr.Dataset:
    """Use the first native cells even when explicitly selected snapshots start later."""
    if source.start in cells.time.values:
        return cells.sel(time=source.start)
    with open_processed_cells(source.run, data_root) as native:
        return native.isel(time=0)[[POSITION, "rr_protein", *DIAGNOSTICS]].load()


def _diagnostics(
    sources: list[PlotSource],
    cells: list[xr.Dataset],
    directory: Path,
    options: PlotOptions,
) -> None:
    chosen: dict[int | str, int] = {}
    for index, source in enumerate(sources):
        key = (
            source.run.condition.pattern_key
            if source.run.condition
            else source.run.run_id
        )
        if key not in chosen or (
            source.dose is not None and source.dose > (sources[chosen[key]].dose or 0)
        ):
            chosen[key] = index
    indices = list(chosen.values())
    fields = ("rr_protein", *DIAGNOSTICS)
    art = figure((6 * 2.7, len(indices) * 2.5))
    axes = art.subplots(len(indices), 6, squeeze=False)
    maximum_volume = 0.0
    try:
        for row, index in enumerate(indices):
            source, data = sources[index], cells[index]
            _, color = _style(source, index, len(sources))
            position = data[POSITION].values
            axes[row, 0].hist(
                position, bins=10, color=color, edgecolor="black", alpha=0.7
            )
            axes[row, 0].set_ylabel("Frequency [-]", **LABEL_STYLE)
            for column, field in enumerate(fields, start=1):
                values = data[field].values * DIAGNOSTIC_FACTORS.get(field, 1.0)
                axes[row, column].plot(
                    position,
                    values,
                    color=color,
                    marker="o",
                    linestyle="",
                    markeredgecolor="black",
                    alpha=0.7,
                )
                axes[row, column].set_ylabel(LABELS[field], **LABEL_STYLE)
                if field in DIAGNOSTIC_FACTORS:
                    maximum_volume = max(maximum_volume, float(np.nanmax(values)))
        for row in range(len(indices)):
            for column in range(6):
                axis = axes[row, column]
                axis.set_xlim(0, 1)
                axis.set_xticks([0, 0.5, 1], labels=["PP", "0.5", "PV"])
                axis.tick_params(axis="x", labelbottom=row == len(indices) - 1)
                if row == len(indices) - 1:
                    axis.set_xlabel("Position [-]", **LABEL_STYLE)
                if column < 5:
                    axis.set_ylim(bottom=0)
                if column in (3, 4):
                    upper = maximum_volume or 1.0
                    axis.set_ylim(-0.05 * upper, 1.05 * upper)
        save_figure(
            art,
            directory / "initial-properties",
            RenderOptions(formats=options.formats, dpi=options.dpi),
        )
    finally:
        art.clear()


def _endpoints(
    aligned: xr.Dataset,
    sources: list[PlotSource],
    directory: Path,
    options: PlotOptions,
    symlog: bool,
    cells: list[xr.Dataset],
) -> None:
    if sources[0].family == "convergence":
        _resolution(
            cells,
            sources,
            directory,
            options,
            symlog,
            tuple(str(field) for field in aligned.field.values),
        )
        return
    fields = [str(field) for field in aligned.field.values]
    art = figure((2.7 * len(fields), 2.5))
    axes = art.subplots(1, len(fields), squeeze=False)[0]
    groups = sorted(
        {source.run.condition.pattern_key for source in sources if source.run.condition}
    ) or [None]
    try:
        for group in groups:
            selected = [
                index
                for index, source in enumerate(sources)
                if group is None
                or (source.run.condition and source.run.condition.pattern_key == group)
            ]
            x = np.array(
                [
                    cast(int, sources[index].dose) + 1
                    if sources[index].dose is not None
                    else int(aligned.num_cells.isel(run=index))
                    for index in selected
                ]
            )
            order = np.argsort(x)
            for axis, field in zip(axes, fields, strict=True):
                y = aligned["mean"].isel(time=-1, run=selected).sel(
                    field=field
                ).values * _factor(field)
                axis.plot(
                    x[order],
                    y[order],
                    marker="o",
                    markersize=6,
                    markeredgecolor="black",
                    color=simulation_color(group, 8) if group is not None else "black",
                    label=pattern_idx2name[group] if group is not None else "Mesh mean",
                )
                axis.set_ylabel(LABELS[field], **LABEL_STYLE)
                axis.set_xlabel(
                    "Dose index" if group is not None else "# Elements", **LABEL_STYLE
                )
                if group is None:
                    axis.set_xscale("log")
        for axis, field in zip(axes, fields, strict=True):
            _scale(axis, field, symlog)
        handles, labels = axes[0].get_legend_handles_labels()
        art.legend(
            handles,
            labels,
            loc="outside lower center",
            ncol=min(6, len(labels)),
            fontsize=8,
        )
        art.suptitle(
            f"Matched-time means at {float(aligned.time.values[-1]) / 3600:g} h"
            + (" | symlog concentration axes" if symlog else "")
        )
        save_figure(
            art,
            directory / ("endpoints-symlog" if symlog else "endpoints"),
            RenderOptions(formats=options.formats, dpi=options.dpi),
        )
    finally:
        art.clear()


def _resolution(
    cells: list[xr.Dataset],
    sources: list[PlotSource],
    directory: Path,
    options: PlotOptions,
    symlog: bool,
    fields: tuple[str, ...],
) -> None:
    """Compare cell means at the historical 1.5 h snapshot when available."""
    times = cells[0].time.values
    time = 5400.0 if 5400.0 in times else float(times[len(times) // 2])
    art = figure((2.7 * len(fields), 2.5))
    axes = art.subplots(1, len(fields), squeeze=False)[0]
    try:
        sizes = np.array([data.sizes["cell"] for data in cells])
        order = np.argsort(sizes)
        for axis, field in zip(axes, fields, strict=True):
            means = []
            for index, (source, data) in enumerate(zip(sources, cells, strict=True)):
                selected = data.sel(time=time)
                mean, _ = weighted_statistics(
                    selected[field].values[None, :], selected[VOLUME].values
                )
                value = float(mean[0]) * _factor(field)
                means.append(value)
                _, color = _style(source, index, len(sources))
                axis.plot(
                    sizes[index],
                    value,
                    marker="o",
                    markersize=6,
                    markeredgecolor="black",
                    color=color,
                    linestyle="",
                )
            axis.plot(sizes[order], np.asarray(means)[order], color="black")
            axis.set_xlabel("# Elements", **LABEL_STYLE)
            axis.set_ylabel(LABELS[field], **LABEL_STYLE)
            axis.set_xlim(left=0, right=axis.get_xlim()[1] * 1.05)
            axis.set_ylim(bottom=0, top=max(axis.get_ylim()[1] * 1.05, 1e-12))
            if symlog:
                _scale(axis, field, True)
        save_figure(
            art,
            directory / ("endpoints-symlog" if symlog else "endpoints"),
            RenderOptions(formats=options.formats, dpi=options.dpi),
        )
    finally:
        art.clear()


def _errors(
    means: pd.DataFrame,
    spatial: pd.DataFrame,
    aligned: xr.Dataset,
    directory: Path,
    options: PlotOptions,
    symlog: bool,
) -> None:
    fields = [str(field) for field in aligned.field.values]
    art = figure((max(8.0, 3.1 * len(fields)), 6.0))
    axes = art.subplots(2, len(fields), squeeze=False)
    try:
        sizes = {
            str(run): int(aligned.num_cells.sel(run=run)) for run in aligned.run.values
        }
        for row, (frame, column, label) in enumerate(
            (
                (means, "absolute_difference", "Mean absolute difference"),
                (spatial, "profile_rmse", "Profile RMS difference"),
            )
        ):
            final = frame[frame.time == frame.time.max()]
            for index, field in enumerate(fields):
                selected = final[final.field == field].copy()
                selected["cells"] = selected.run.map(sizes)
                selected = selected.sort_values("cells")
                axes[row, index].plot(
                    selected.cells,
                    selected[column] * _factor(field),
                    marker="o",
                    markersize=3,
                )
                axes[row, index].set_xscale("log")
                axes[row, index].set_xlabel("Lobulus cells")
                _scale(axes[row, index], "difference", symlog)
                if row == 0:
                    axes[row, index].set_title(LABELS[field], fontsize=10)
            axes[row, 0].set_ylabel(label)
        art.suptitle(
            f"Difference from finest available mesh ({means.reference_run.iloc[0]}); matched times and bins"
        )
        save_figure(
            art,
            directory
            / ("convergence-errors-symlog" if symlog else "convergence-errors"),
            RenderOptions(formats=options.formats, dpi=options.dpi),
        )
    finally:
        art.clear()


def plot_configuration(
    sources: list[PlotSource], options: PlotOptions
) -> dict[str, object]:
    """Describe plot choices independently of numerical reduction options."""
    return {
        "version": ANALYSIS_PLOT_VERSION,
        "study": sources[0].run.study,
        "runs": [source.run.run_id for source in sources],
        "formats": list(options.formats),
        "dpi": options.dpi,
        "spread": options.spread,
        "symlog": options.convergence_symlog,
    }


def cached_analysis_figures(
    sources: list[PlotSource], output_root: Path, options: PlotOptions
) -> bool:
    """Check all study exports without opening numerical data or figures."""
    directory = sources[0].directory(output_root, "analysis").parent
    return (
        current_manifest(directory, plot_configuration(sources, options), output_root)
        is not None
    )


def render_analysis(
    sources: list[PlotSource],
    prepared: list[dict[str, object]],
    output_root: Path,
    options: PlotOptions,
    *,
    force: bool = False,
) -> Path:
    """Publish self-contained tables, NetCDF comparisons and their figures atomically."""
    options.validate()
    directory = sources[0].directory(output_root, "analysis").parent
    paths = tuple(
        source.directory(output_root, "analysis-data") / "manifest.json"
        for source in sources
    )
    configuration = plot_configuration(sources, options)
    with run_lock(directory):
        cached = current_manifest(directory, configuration, output_root)
        if cached and not force:
            return generation_path(directory, cached)
        before = [fingerprint(path) for path in paths]
        native = [load_analysis(source, output_root, "native") for source in sources]
        profiles = [
            load_analysis(source, output_root, "profiles") for source in sources
        ]
        cells = [load_analysis(source, output_root, "cells") for source in sources]
        initial = [
            _initial_cells(
                source,
                data,
                Path(str(cast(dict[str, object], item["metadata"])["data_root"])),
            )
            for source, data, item in zip(sources, cells, prepared, strict=True)
        ]
        aligned_parts = [
            load_analysis(source, output_root, "aligned") for source in sources
        ]
        run_ids = [source.run.run_id for source in sources]
        aligned = xr.concat(aligned_parts, dim=xr.IndexVariable("run", run_ids))
        aligned = aligned.assign_coords(
            num_cells=(
                "run",
                [
                    cast(int, cast(dict[str, object], item["metadata"])["num_cells"])
                    for item in prepared
                ],
            )
        )
        all_profiles = xr.concat(profiles, dim=xr.IndexVariable("run", run_ids))
        with tempfile.TemporaryDirectory(
            prefix=".partial-", dir=directory
        ) as temporary:
            staging = Path(temporary)
            tables = []
            for source, dataset in zip(sources, native, strict=True):
                frame = dataset.to_dataframe().reset_index()
                frame.insert(0, "run", source.run.run_id)
                frame["study"] = source.run.study
                frame["pattern"] = (
                    source.run.condition.pattern_key if source.run.condition else None
                )
                frame["dose"] = source.dose + 1 if source.dose is not None else None
                tables.append(frame)
            pd.concat(tables, ignore_index=True).to_csv(
                staging / "native-timecourses.tsv", sep="\t", index=False
            )
            aligned.to_dataframe().reset_index().to_csv(
                staging / "matched-timecourses.tsv", sep="\t", index=False
            )
            all_profiles.to_dataframe().reset_index().to_csv(
                staging / "spatial-profiles.tsv", sep="\t", index=False
            )
            result = xr.merge(
                [
                    aligned.rename(
                        {
                            "field": "analysis_field",
                            **{name: f"aligned_{name}" for name in aligned.data_vars},
                        }
                    ),
                    all_profiles.rename(
                        {
                            "time": "spatial_time",
                            "position": "position_bin",
                            **{
                                name: f"profile_{name}"
                                for name in all_profiles.data_vars
                            },
                        }
                    ),
                ]
            )
            result = result.assign_coords(
                native_end_s=("run", [source.end for source in sources]),
                pattern=(
                    "run",
                    [
                        source.run.condition.pattern_key if source.run.condition else -1
                        for source in sources
                    ],
                ),
                dose=(
                    "run",
                    [
                        source.dose + 1 if source.dose is not None else -1
                        for source in sources
                    ],
                ),
            )
            result.attrs.update(
                definition="Cell-volume-weighted population mean/std; source scalar units; native histories in native-timecourses.tsv",
                spatial_time_units="s",
                comparison="matched continuous cell interpolation with previous discrete values; no extrapolation",
            )
            result.to_netcdf(staging / "analysis.nc", engine="scipy")
            means, spatial_errors = (
                comparisons(aligned, all_profiles, sources)
                if sources[0].family == "convergence"
                else (None, None)
            )
            if means is not None and spatial_errors is not None:
                means.to_csv(staging / "convergence-means.tsv", sep="\t", index=False)
                spatial_errors.to_csv(
                    staging / "convergence-profiles.tsv", sep="\t", index=False
                )
            with mpl.rc_context({"svg.fonttype": "none", "pdf.fonttype": 42}):
                for symlog in (
                    (False, True)
                    if sources[0].family == "convergence" and options.convergence_symlog
                    else (False,)
                ):
                    suffix = "-symlog" if symlog else ""
                    _curves(
                        sources,
                        aligned_parts,
                        staging,
                        f"timecourses{suffix}",
                        options,
                        tuple(str(field) for field in aligned.field.values),
                        symlog=symlog,
                    )
                    for index, time in enumerate(profiles[0].time.values):
                        _curves(
                            sources,
                            cells,
                            staging,
                            f"spatial-{index:02}{suffix}",
                            options,
                            tuple(str(field) for field in aligned.field.values),
                            spatial_time=float(time),
                            symlog=symlog,
                        )
                    _endpoints(aligned, sources, staging, options, symlog, cells)
                    if means is not None and spatial_errors is not None:
                        _errors(
                            means, spatial_errors, aligned, staging, options, symlog
                        )
                _diagnostics(sources, initial, staging, options)
            atomic_json(
                staging / "provenance.json",
                {
                    **configuration,
                    "analyses": prepared,
                    "sources": [fingerprint(path) for path in paths],
                },
            )
            links = [
                f'<li><a href="{html.escape(path.name, quote=True)}">{html.escape(path.name)}</a></li>'
                for path in sorted(staging.iterdir())
                if path.suffix in {".nc", ".tsv", ".json"}
            ]
            cards = []
            for name in sorted(
                {
                    path.stem
                    for path in staging.iterdir()
                    if path.suffix in {".png", ".svg", ".pdf"}
                }
            ):
                exports = " | ".join(
                    f'<a href="{name}.{extension}">{extension.upper()}</a>'
                    for extension in options.formats
                )
                preview = (
                    f'<img loading="lazy" src="{name}.png" style="width:100%">'
                    if "png" in options.formats
                    else ""
                )
                cards.append(
                    f"<article><h2>{name}</h2>{preview}<p>{exports}</p></article>"
                )
            atomic_text(
                staging / "index.html",
                f'<!doctype html><html lang="en"><meta charset="utf-8"><title>SPT analysis</title><style>body{{font:16px system-ui;margin:2rem}}article{{max-width:1500px;margin:2rem 0}}h2{{font-size:1.1rem}}</style><h1>{sources[0].run.study}</h1><p>Volume-weighted means; SD describes spatial variation. Concentrations retain source units. Necrosis figures use percent; serialized values remain fractions. Empty spatial bins remain missing.</p><ul>{"".join(links)}</ul>{"".join(cards)}</html>',
            )
            generation = publish(
                directory,
                staging,
                configuration,
                paths,
                before,
                {"runs": len(sources), "figures": len(cards)},
            )
        for dataset in (
            *native,
            *profiles,
            *initial,
            *cells,
            *aligned_parts,
            aligned,
            all_profiles,
            result,
        ):
            dataset.close()
        return generation
