"""Compose SPT time, position, dose and mesh comparisons from serialized results."""

from __future__ import annotations

import html
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import cast

import matplotlib as mpl
import numpy as np
import pandas as pd
import xarray as xr
from matplotlib.axes import Axes
from matplotlib.colors import to_hex
from matplotlib.ticker import FixedLocator, SymmetricalLogLocator

from research.spt.analysis_data import DIAGNOSTICS, VOLUME, load_analysis
from research.spt.artifacts import current_manifest, publish
from research.spt.cache import _generation, _run_lock, fingerprint
from research.spt.cli import atomic_json
from research.spt.conditions import pattern_idx2name
from research.spt.plot_data import PlotSource
from research.spt.plot_render import LAYERS, RenderOptions, figure, save_figure

ANALYSIS_PLOT_VERSION = 3
LABELS = {layer.sid: layer.title for layer in LAYERS}
LABELS["rr_necrosis"] = "Necrotic volume [%]"
LABELS.update(
    {
        VOLUME: "Cell volume [source units]",
        "rr_Vext": "Fluid volume [source units]",
        "rr_Vli": "Solid volume [source units]",
        "pressure": "Pressure [source units]",
    }
)


@dataclass(frozen=True)
class PlotOptions:
    """Select compact figures and optional population spread or symlog comparisons."""

    formats: tuple[str, ...] = ("svg", "pdf", "png")
    dpi: int = 110
    spread: bool = False
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
        dose = source.run.condition.boundary_flow_key + 1
        return f"Dose {dose}", to_hex(
            mpl.colormaps["viridis"](0.08 + 0.85 * (dose - 1) / 7)
        )
    return source.run.run_id.removeprefix("convergence_"), to_hex(
        mpl.colormaps["viridis"](index / max(count - 1, 1))
    )


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
    axis.grid(alpha=0.18)


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
) -> None:
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
    art = figure((max(8.0, 3.0 * len(fields)), 2.35 * len(patterns) + 1.2))
    axes = art.subplots(
        len(patterns), len(fields), squeeze=False, sharex=True, sharey="col"
    )
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
                if spatial_time is None:
                    x, y, std = (
                        dataset.time.values / 3600,
                        dataset[f"{field}_mean"].values,
                        dataset[f"{field}_std"].values,
                    )
                else:
                    values = dataset.sel(time=spatial_time, field=field)
                    x, y, std = (
                        dataset.position.values,
                        values["mean"].values,
                        values["std"].values,
                    )
                factor = _factor(field)
                (line,) = axis.plot(
                    x, y * factor, color=color, linewidth=1.15, label=label
                )
                if options.spread:
                    axis.fill_between(
                        x,
                        (y - std) * factor,
                        (y + std) * factor,
                        color=color,
                        alpha=0.10,
                    )
                legend[label] = line
        for row in range(len(patterns)):
            for column, field in enumerate(fields):
                axis = axes[row, column]
                _scale(axis, field, symlog)
                if row == 0:
                    axis.set_title(LABELS[field], fontsize=10)
                if row == len(patterns) - 1:
                    axis.set_xlabel(
                        "Time [h]"
                        if spatial_time is None
                        else "Normalized position (PP to PV)"
                    )
                if spatial_time is not None:
                    axis.set_xlim(0, 1)
            if simulation:
                axes[row, 0].set_ylabel(pattern_idx2name[patterns[row]], fontsize=10)
        art.legend(
            list(legend.values()),
            list(legend),
            loc="outside lower center",
            ncol=min(8, len(legend)),
            fontsize=9,
        )
        suffix = " | symlog concentration axes" if symlog else ""
        time_label = (
            "Native-time volume-weighted means"
            if spatial_time is None
            else f"Volume-weighted spatial profiles at {spatial_time / 3600:g} h"
        )
        spread_note = (
            "\nShading: volume-weighted population SD across cells; not a confidence interval"
            if options.spread
            else ""
        )
        art.suptitle(
            f"{sources[0].family.title()} | {sources[0].geometry} | {time_label}{suffix}{spread_note}",
            fontsize=12,
        )
        save_figure(
            art,
            directory / name,
            RenderOptions(formats=options.formats, dpi=options.dpi),
        )
    finally:
        art.clear()


def _diagnostics(
    sources: list[PlotSource],
    profiles: list[xr.Dataset],
    directory: Path,
    options: PlotOptions,
) -> None:
    chosen = {}
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
    art = figure((17, 2.15 * len(indices) + 0.8))
    axes = art.subplots(len(indices), 6, squeeze=False, sharex=True, sharey="col")
    try:
        for row, index in enumerate(indices):
            source, data = sources[index], profiles[index].isel(time=0)
            volume = data.volume.values
            axes[row, 0].bar(
                data.position.values,
                volume / volume.sum(),
                width=0.9 / len(volume),
                color="tab:blue",
            )
            axes[row, 0].set_ylabel(
                pattern_idx2name[source.run.condition.pattern_key]
                if source.run.condition
                else source.run.run_id,
                fontsize=9,
            )
            for column, field in enumerate(fields, start=1):
                values = data.sel(field=field)
                axes[row, column].plot(data.position, values["mean"], color="tab:blue")
            for column in range(6):
                axis = axes[row, column]
                axis.grid(alpha=0.18)
                axis.set_xlim(0, 1)
                if row == 0:
                    axis.set_title(
                        "Volume fraction per bin"
                        if column == 0
                        else LABELS[fields[column - 1]],
                        fontsize=10,
                    )
                if row == len(indices) - 1:
                    axis.set_xlabel("Normalized position")
        art.suptitle(
            f"Initial properties | {sources[0].run.study} | highest selected dose per pattern",
            fontsize=12,
        )
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
) -> None:
    fields = [str(field) for field in aligned.field.values]
    art = figure((max(8.0, 3.1 * len(fields)), 3.6))
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
                    markersize=3,
                    label=pattern_idx2name[group] if group is not None else "Mesh mean",
                )
                axis.set_title(LABELS[field], fontsize=10)
                axis.set_xlabel("Dose index" if group is not None else "Lobulus cells")
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
    with _run_lock(directory):
        cached = current_manifest(directory, configuration, output_root)
        if cached and not force:
            return _generation(directory, cached)
        before = [fingerprint(path) for path in paths]
        native = [load_analysis(source, output_root, "native") for source in sources]
        profiles = [
            load_analysis(source, output_root, "profiles") for source in sources
        ]
        initial = [load_analysis(source, output_root, "initial") for source in sources]
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
                        native,
                        staging,
                        f"timecourses{suffix}",
                        options,
                        tuple(str(field) for field in aligned.field.values),
                        symlog=symlog,
                    )
                    for index, time in enumerate(profiles[0].time.values):
                        _curves(
                            sources,
                            profiles,
                            staging,
                            f"spatial-{index:02}{suffix}",
                            options,
                            tuple(str(field) for field in aligned.field.values),
                            spatial_time=float(time),
                            symlog=symlog,
                        )
                    _endpoints(aligned, sources, staging, options, symlog)
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
            (staging / "index.html").write_text(
                f'<!doctype html><html lang="en"><meta charset="utf-8"><title>SPT analysis</title><style>body{{font:16px system-ui;margin:2rem}}article{{max-width:1500px;margin:2rem 0}}h2{{font-size:1.1rem}}</style><h1>{sources[0].run.study}</h1><p>Volume-weighted means; SD describes spatial variation. Concentrations retain source units. Necrosis figures use percent; serialized values remain fractions. Empty spatial bins remain missing.</p><ul>{"".join(links)}</ul>{"".join(cards)}</html>',
                encoding="utf-8",
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
            *aligned_parts,
            aligned,
            all_profiles,
            result,
        ):
            dataset.close()
        return generation
