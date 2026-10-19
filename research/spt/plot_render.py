"""Render reusable planar figures and stream requested videos without PNG panels."""

from __future__ import annotations

import json
import shutil
import tempfile
import time
from collections.abc import Sequence
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Literal, cast

import matplotlib as mpl
import numpy as np
from matplotlib.animation import FFMpegWriter
from matplotlib.axes import Axes
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.collections import PolyCollection
from matplotlib.colors import SymLogNorm
from matplotlib.figure import Figure
from matplotlib.ticker import FixedLocator
from numpy.typing import ArrayLike

from research.spt.cache import (
    ProcessingResult,
    _generation,
    _run_lock,
    _unchanged,
    fingerprint,
)
from research.spt.catalog import CatalogValue
from research.spt.cli import atomic_json
from research.spt.plot_data import PlotSource, current_manifest, publish

RENDER_VERSION = 3


@dataclass(frozen=True)
class Layer:
    """Describe scalar labels, colors and optional explicit display limits."""

    sid: str
    title: str
    colormap: str = "magma"
    limits: tuple[float, float] | None = None
    factor: float = 1.0
    scale: Literal["linear", "symlog"] = "linear"


LAYERS = (
    Layer("rr_(S_ext)", "Substrate [mM]"),
    Layer("rr_(P_ext)", "Product [mM]", "viridis"),
    Layer("rr_protein", "Protein", "viridis"),
    Layer("rr_(T)", "Toxicant [mM]"),
    Layer("rr_necrosis", "Necrosis (fraction)", "Reds", (0.0, 1.0)),
)


@dataclass(frozen=True)
class RenderOptions:
    """Select outputs explicitly; videos and standalone panels are opt-in."""

    layers: tuple[Layer, ...] = LAYERS
    formats: tuple[str, ...] = ("svg", "pdf", "png")
    snapshots: tuple[float, ...] = ()
    video_times: tuple[float, ...] = ()
    dpi: int = 110
    fps: int = 12
    panels: bool = False
    gif: bool = False

    def validate(self) -> None:
        """Reject invalid rendering choices before creating any output."""
        if not self.layers or len({layer.sid for layer in self.layers}) != len(
            self.layers
        ):
            raise ValueError("Select nonempty unique scalar layers")
        if not self.formats or set(self.formats) - {"svg", "pdf", "png"}:
            raise ValueError("Figure formats must be svg, pdf or png")
        if self.dpi < 1 or self.fps < 1:
            raise ValueError("dpi and fps must be positive")
        if self.gif and not self.video_times:
            raise ValueError("GIF output requires video frames")
        if self.video_times and shutil.which("ffmpeg") is None:
            raise FileNotFoundError("FFmpeg is required only when videos are requested")
        for layer in self.layers:
            if layer.scale not in {"linear", "symlog"}:
                raise ValueError(f"Unknown color scale: {layer.scale}")
            if not np.isfinite(layer.factor) or (
                layer.limits is not None
                and (
                    not np.isfinite(layer.limits).all()
                    or layer.limits[0] >= layer.limits[1]
                )
            ):
                raise ValueError(f"Invalid display scale: {layer.sid}")
            mpl.colormaps[layer.colormap]


@dataclass
class Samples:
    """Own compact selected arrays after their NPZ file has been closed."""

    faces: np.ndarray
    cell_ids: np.ndarray
    values: np.ndarray
    sectors: np.ndarray
    times: np.ndarray

    @classmethod
    def load(cls, path: Path) -> Samples:
        """Load numerical arrays without allowing pickle payloads."""
        with np.load(path, allow_pickle=False) as data:
            return cls(
                *(
                    np.asarray(data[name])
                    for name in ("faces", "cell_ids", "values", "sectors", "times")
                )
            )

    def frame(self, time_value: float) -> np.ndarray:
        """Select an already sampled frame exactly, without a second interpolation."""
        index = int(np.searchsorted(self.times, time_value))
        if index >= len(self.times) or self.times[index] != time_value:
            raise ValueError(f"Frame was not prepared: {time_value:g} s")
        return self.values[index]


def display_limits(
    layer: Layer, ranges: dict[str, tuple[float, float]]
) -> tuple[float, float]:
    """Apply consistent conversions and pad constant ranges for valid normalization."""
    if layer.limits is not None:
        return layer.limits
    low, high = sorted(value * layer.factor for value in ranges[layer.sid])
    if low == high:
        if low > 0:
            low = 0.0
        elif high < 0:
            high = 0.0
        else:
            high = 1.0
    return low, high


def figure(size: tuple[float, float]) -> Figure:
    """Create a headless figure independent of pyplot and a GUI session."""
    result = Figure(figsize=size, layout="constrained")
    FigureCanvasAgg(result)
    return result


def panel(
    axis: Axes,
    samples: Samples,
    layer: Layer,
    values: np.ndarray,
    limits: tuple[float, float],
    *,
    mesh: bool = False,
) -> PolyCollection:
    """Draw the original cell values on planar faces without point averaging."""
    collection = PolyCollection(
        cast(Sequence[ArrayLike], samples.faces),
        cmap=layer.colormap,
        linewidths=min(0.13, 5.0 / np.sqrt(len(samples.faces))) if mesh else 0,
        edgecolors="black" if mesh else "none",
        antialiaseds=False,
    )
    collection.set_array(values[samples.cell_ids] * layer.factor)
    if layer.scale == "symlog":
        collection.set_norm(
            SymLogNorm(
                linthresh=max(abs(value) for value in limits) * 1e-6,
                vmin=limits[0],
                vmax=limits[1],
            )
        )
    else:
        collection.set_clim(*limits)
    axis.add_collection(collection)
    low = samples.faces.min(axis=(0, 1))
    high = samples.faces.max(axis=(0, 1))
    center = (low + high) / 2
    extent = max(float((high - low).max()), np.finfo(float).eps) * 0.56
    axis.set_xlim(float(center[0] - extent), float(center[0] + extent))
    axis.set_ylim(float(center[1] - extent), float(center[1] + extent))
    axis.set_aspect("equal")
    axis.set_axis_off()
    return collection


def label_sectors(axis: Axes, samples: Samples) -> None:
    """Number independent sectors once, without covering every time-course panel."""
    centers = samples.faces.mean(axis=1)
    sectors = samples.sectors[samples.cell_ids]
    for sector in np.unique(sectors):
        if sector < 0:
            continue
        position = centers[sectors == sector].mean(axis=0)
        axis.text(
            float(position[0]),
            float(position[1]),
            str(int(sector) + 1),
            ha="center",
            va="center",
            fontsize=8,
            bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "none", "pad": 1},
        )


def save_figure(art: Figure, destination: Path, options: RenderOptions) -> None:
    """Save true vector SVG/PDF and optional previews from the same figure."""
    # Matplotlib otherwise rasterizes long colorbar gradients in vector exports.
    for axis in art.axes:
        for collection in axis.collections:
            collection.set_rasterized(False)
            norm = collection.norm
            if isinstance(norm, SymLogNorm):
                for dimension in (axis.xaxis, axis.yaxis):
                    if dimension.get_scale() == "symlog":
                        low, high = dimension.get_view_interval()
                        dimension.set_major_locator(
                            FixedLocator(
                                [
                                    value
                                    for value in dimension.get_majorticklocs()
                                    if low <= value <= high
                                    and (value == 0 or abs(value) >= norm.linthresh)
                                ]
                            )
                        )
    destination.parent.mkdir(parents=True, exist_ok=True)
    for extension in options.formats:
        art.savefig(destination.with_suffix(f".{extension}"), dpi=options.dpi)


def render_timecourse(
    source: PlotSource,
    samples: Samples,
    directory: Path,
    options: RenderOptions,
    ranges: dict[str, tuple[float, float]],
) -> None:
    """Render exactly the requested times and fields in one composed figure."""
    art = figure((3.0 * len(options.layers), 2.6 * len(options.snapshots) + 0.9))
    axes = art.subplots(len(options.snapshots), len(options.layers), squeeze=False)
    try:
        for row, time_value in enumerate(options.snapshots):
            values = samples.frame(time_value)
            for column, layer in enumerate(options.layers):
                axis = axes[row, column]
                collection = panel(
                    axis, samples, layer, values[column], display_limits(layer, ranges)
                )
                if row == 0:
                    axis.set_title(layer.title, fontsize=11)
                if column == 0:
                    axis.text(
                        -0.03,
                        0.5,
                        f"{time_value / 3600:g} h",
                        transform=axis.transAxes,
                        ha="right",
                        va="center",
                        rotation=90,
                        fontsize=10,
                    )
                if row == 0 and column == 0 and source.sector_labels:
                    label_sectors(axis, samples)
                if row == len(options.snapshots) - 1:
                    art.colorbar(
                        collection,
                        ax=list(axes[:, column]),
                        location="bottom",
                        shrink=0.9,
                        aspect=20,
                        pad=0.01,
                    )
        art.suptitle(source.title, fontsize=13)
        if source.sector_labels:
            art.supxlabel(
                " | ".join(
                    f"{index + 1}: {name}"
                    for index, name in enumerate(source.sector_labels)
                ),
                fontsize=9,
            )
        save_figure(art, directory / "timecourse", options)
    finally:
        art.clear()


def render_panels(
    source: PlotSource,
    samples: Samples,
    directory: Path,
    options: RenderOptions,
    ranges: dict[str, tuple[float, float]],
) -> None:
    """Produce individual field panels only when they have been requested."""
    for index, layer in enumerate(options.layers):
        art = figure((4.0, 4.0))
        try:
            axis = art.subplots()
            collection = panel(
                axis,
                samples,
                layer,
                samples.frame(options.snapshots[0])[index],
                display_limits(layer, ranges),
            )
            art.colorbar(collection, ax=axis, shrink=0.75)
            for frame, time_value in enumerate(options.snapshots):
                collection.set_array(
                    samples.frame(time_value)[index][samples.cell_ids] * layer.factor
                )
                axis.set_title(f"{layer.title} | {time_value / 3600:g} h", fontsize=11)
                save_figure(
                    art,
                    directory / "panels" / f"field_{index + 1:02}_t{frame:05}",
                    options,
                )
        finally:
            art.clear()


def render_video(
    samples: Samples,
    directory: Path,
    options: RenderOptions,
    ranges: dict[str, tuple[float, float]],
) -> None:
    """Reuse artists and stream a single combined frame directly into FFmpeg."""
    art = figure((3.2 * len(options.layers), 3.8))
    axes = art.subplots(1, len(options.layers), squeeze=False)[0]
    try:
        collections = []
        first = samples.frame(options.video_times[0])
        for index, (axis, layer) in enumerate(zip(axes, options.layers, strict=True)):
            collection = panel(
                axis, samples, layer, first[index], display_limits(layer, ranges)
            )
            axis.set_title(layer.title, fontsize=10)
            art.colorbar(collection, ax=axis, shrink=0.75)
            collections.append(collection)
        title = art.suptitle("")
        writer = FFMpegWriter(
            fps=options.fps,
            codec="h264",
            extra_args=["-pix_fmt", "yuv420p", "-movflags", "+faststart"],
        )
        with writer.saving(art, str(directory / "timecourse.mp4"), options.dpi):
            for time_value in options.video_times:
                values = samples.frame(time_value)
                for index, (collection, layer) in enumerate(
                    zip(collections, options.layers, strict=True)
                ):
                    collection.set_array(values[index][samples.cell_ids] * layer.factor)
                title.set_text(f"{time_value / 3600:g} h")
                writer.grab_frame(facecolor="white")
        if options.gif:
            from porous_media.visualization.video import create_gif_from_video

            create_gif_from_video(
                directory / "timecourse.mp4",
                directory / "timecourse.gif",
                overwrite=False,
            )
    finally:
        art.clear()


def render_source(
    source: PlotSource,
    prepared: dict[str, object],
    output_root: Path,
    options: RenderOptions,
    ranges: dict[str, tuple[float, float]],
    *,
    force: bool = False,
    dry_run: bool = False,
) -> ProcessingResult:
    """Publish one complete figure set, preserving a previous generation on failure."""
    started = time.perf_counter()
    try:
        output_root = output_root.expanduser().resolve()
        options.validate()
        if not _unchanged(prepared.get("sources"), source.run.root):
            raise ValueError("Plot samples are stale; prepare them again")
        sample_directory = source.directory(output_root, "plot-data")
        sample_generation = _generation(sample_directory, prepared)
        paths = tuple(
            sorted(path for path in sample_generation.iterdir() if path.is_file())
        )
        directory = source.directory(output_root, "figures")
        configuration: dict[str, object] = {
            "version": RENDER_VERSION,
            "source": str(source.run.xdmf),
            "title": source.title,
            "sector_labels": list(source.sector_labels),
            "options": json.loads(json.dumps(asdict(options))),
            "ranges": {name: list(values) for name, values in ranges.items()},
        }
        cached = current_manifest(directory, configuration, output_root)
        if dry_run:
            return ProcessingResult(
                source.run.study,
                source.run.run_id,
                "cached" if cached and not force else "pending",
                "cached" if cached and not force else "render required",
                0.0,
                {},
            )
        with _run_lock(directory):
            cached = current_manifest(directory, configuration, output_root)
            if cached and not force:
                return ProcessingResult(
                    source.run.study,
                    source.run.run_id,
                    "cached",
                    "cached",
                    time.perf_counter() - started,
                    figure_row(source, _generation(directory, cached), options),
                )
            before = [fingerprint(path) for path in paths]
            with _run_lock(sample_directory):
                if not _unchanged(prepared.get("sources"), source.run.root):
                    raise ValueError("Plot samples are stale; prepare them again")
                samples = Samples.load(sample_generation / "samples.npz")
            with tempfile.TemporaryDirectory(
                prefix=".partial-", dir=directory
            ) as temporary:
                staging = Path(temporary)
                with mpl.rc_context({"svg.fonttype": "none", "pdf.fonttype": 42}):
                    render_timecourse(source, samples, staging, options, ranges)
                    if options.panels:
                        render_panels(source, samples, staging, options, ranges)
                    if options.video_times:
                        render_video(samples, staging, options, ranges)
                atomic_json(
                    staging / "provenance.json",
                    {**configuration, "sampling": prepared["metadata"]},
                )
                if before != [fingerprint(path) for path in paths] or not _unchanged(
                    prepared.get("sources"), source.run.root
                ):
                    raise RuntimeError("Plot inputs changed during rendering")
                metadata: dict[str, object] = {
                    "study": source.run.study,
                    "run_id": source.run.run_id,
                    "geometry": source.geometry,
                    "num_snapshots": len(options.snapshots),
                    "num_fields": len(options.layers),
                }
                generation = publish(
                    directory, staging, configuration, paths, before, metadata
                )
                row = figure_row(source, generation, options)
            return ProcessingResult(
                source.run.study,
                source.run.run_id,
                "rendered",
                "rendered",
                time.perf_counter() - started,
                row,
            )
    except Exception as error:
        return ProcessingResult(
            source.run.study,
            source.run.run_id,
            "failed",
            f"{type(error).__name__}: {error}",
            time.perf_counter() - started,
            {},
        )


def figure_row(
    source: PlotSource, generation: Path, options: RenderOptions
) -> dict[str, CatalogValue]:
    """Expose published artifact paths without a second manifest transaction."""
    row: dict[str, CatalogValue] = {
        "study": source.run.study,
        "run_id": source.run.run_id,
        "geometry": source.geometry,
        "num_snapshots": len(options.snapshots),
        "num_fields": len(options.layers),
    }
    for extension in options.formats:
        row[extension] = str(generation / f"timecourse.{extension}")
    if options.video_times:
        row["video"] = str(generation / "timecourse.mp4")
    if options.gif:
        row["gif"] = str(generation / "timecourse.gif")
    return row
