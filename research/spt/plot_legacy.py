"""Keep historical planar plotting calls on the shared selective renderer."""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path
from typing import TYPE_CHECKING

from research.spt.cache import ProcessingResult, validate_roots
from research.spt.catalog import RunSource
from research.spt.plot_data import PlotSource, _time, source_row
from research.spt.plot_geometry import process_plots
from research.spt.plot_render import Layer, RenderOptions

if TYPE_CHECKING:
    from porous_media.visualization.pyvista_visualization import DataLayer


def _sources(paths: Iterable[Path], output_root: Path) -> list[PlotSource]:
    """Materialize an iterable once and reject empty or ambiguous file selections."""
    selected = tuple(Path(path).expanduser().resolve() for path in paths)
    if not selected:
        raise ValueError("Select at least one XDMF history")
    if len({path.stem for path in selected}) != len(selected):
        raise ValueError("Selected descriptors must have unique run names")
    result = []
    for path in selected:
        validate_roots(path.parent, output_root)
        run = RunSource(
            path.parent,
            "simulations_sixth",
            path.stem,
            None,
            "explicit descriptor",
            path,
            path.with_suffix(".feb"),
        )
        row = source_row(run, output_root, output_root)
        result.append(
            PlotSource(
                run, "sixth", _time(row, "tstart_s"), _time(row, "tend_s"), path.stem
            )
        )
    return result


def visualize_spt_2d(
    xdmf_paths: Iterable[Path],
    data_layers: Iterable[DataLayer],
    results_dir: Path,
    selection: Iterable[str],
    create_panels: bool = False,
    *,
    snapshots: int = 5,
    formats: tuple[str, ...] = ("svg", "pdf", "png"),
    video_frames: int = 0,
) -> list[ProcessingResult]:
    """Render selected planar scalar fields without interpolated XDMF copies.

    Historical callers retain their layer labels, colormaps, limits and unit
    conversions. Generator inputs are supported. Standalone panels and videos
    now require explicit requests; composed figures are always generated.
    """
    sources = _sources(xdmf_paths, results_dir)
    chosen = tuple(selection)
    available = {layer.sid: layer for layer in data_layers}
    layers = []
    for name in chosen:
        if name not in available:
            raise ValueError(f"Unknown data layer: {name}")
        layer = available[name]
        if not isinstance(layer.colormap, str):
            raise ValueError("Planar figures require a named Matplotlib colormap")
        layers.append(
            Layer(
                name,
                layer.title,
                layer.colormap,
                layer.color_limits,
                1.0 if layer.conversion_factor is None else layer.conversion_factor,
            )
        )
    results = process_plots(
        sources[0].run.root,
        results_dir,
        results_dir,
        options=RenderOptions(tuple(layers), formats, panels=create_panels),
        snapshots=snapshots,
        video_frames=video_frames,
        workers=1,
        provided_sources=sources,
    )
    failures = [result.message for result in results if result.status == "failed"]
    if failures:
        raise RuntimeError("; ".join(failures))
    return results


def necrosis_plots(
    xdmf_paths: Iterable[Path], results_dir: Path
) -> list[ProcessingResult]:
    """Create the final-state necrosis overview without generating unused panels."""
    sources = _sources(xdmf_paths, results_dir)
    return process_plots(
        sources[0].run.root,
        results_dir,
        results_dir,
        options=RenderOptions(
            (Layer("rr_necrosis", "Necrosis (fraction)", "Reds", (0, 1)),)
        ),
        snapshots=1,
        workers=1,
        per_history=False,
        provided_sources=sources,
    )
