"""Create flat square tissue panels by repeating the linear-increase SPT lobule."""

from __future__ import annotations

import argparse
import json
import sys
import tempfile
from collections.abc import Sequence
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import cast

import matplotlib as mpl
import numpy as np
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.cm import ScalarMappable
from matplotlib.collections import PolyCollection
from matplotlib.colors import Normalize
from matplotlib.figure import Figure
from matplotlib.patches import Patch
from numpy.typing import ArrayLike

from research import RESULTS_DIR
from research.spt.artifacts import (
    atomic_json,
    atomic_text,
    current_manifest,
    fingerprint,
    generation_path,
    publish,
    run_lock,
)
from research.spt.fields import LAYERS, SUMMARY_FIELDS, Layer
from research.spt.plot_render import RenderOptions, display_limits, figure
from research.spt.publication import export_figure, simulation_cohort
from research.spt.publication_data import SIMULATION_TIMES, PublicationData
from research.spt.publication_geometry import save_publication, shared_colorbar
from research.spt.publication_patterns import add_pattern_icon
from research.spt.report import build_report
from research.spt.tissue_geometry import TissuePatch, tile_lobules

TISSUE_VERSION = 2
TISSUE_LAYERS = (
    replace(LAYERS[2], title="CYP2E1"),
    LAYERS[3],
    replace(LAYERS[4], title="Necrosis", limits=(0, 1)),
)
FIELD_NAMES = {
    "rr_protein": "cyp2e1",
    "rr_(T)": "toxic-compound",
    "rr_necrosis": "necrosis",
}


@dataclass(frozen=True)
class TissueOptions:
    """Select a source dose, native times and exact square tissue coverage."""

    dose: int = 6
    times: tuple[float, ...] = SIMULATION_TIMES
    width: float = 2.5
    dpi: int = 300
    formats: tuple[str, ...] = ("svg", "pdf", "png")
    vessel_color: str = "white"
    mesh: bool = False

    def validate(self) -> None:
        """Reject unsupported choices before writing any output."""
        RenderOptions(formats=self.formats, dpi=self.dpi).validate()
        if not 1 <= self.dose <= 8:
            raise ValueError("Dose index must be between 1 and 8")
        if not np.isfinite(self.width) or self.width < 2.5:
            raise ValueError("Square width must be finite and at least 2.5 lobules")
        if (
            not self.times
            or not np.isfinite(self.times).all()
            or min(self.times) < 0
            or np.any(np.diff(self.times) <= 0)
        ):
            raise ValueError(
                "Snapshot times must be finite, nonnegative and increasing"
            )
        if self.vessel_color not in {"white", "lightgray"}:
            raise ValueError("Vessel color must be white or lightgray")


def patch_image(
    patch: TissuePatch,
    values: np.ndarray,
    layer: Layer,
    limits: tuple[float, float],
    options: TissueOptions,
) -> np.ndarray:
    """Render native XY cells directly from above into a square, margin-free tile.

    This deliberately does not use the study's tilted geometry renderer. Cells
    retain their owning scalar value at clipped borders and across every copy.
    The background represents native vessel lumens, never fabricated tissue data.
    """
    art = Figure(figsize=(3, 3), dpi=options.dpi, facecolor=options.vessel_color)
    canvas = FigureCanvasAgg(art)
    try:
        axis = art.add_axes((0, 0, 1, 1), facecolor=options.vessel_color)
        collection = PolyCollection(
            cast(Sequence[ArrayLike], patch.faces),
            cmap=layer.colormap,
            norm=Normalize(*limits),
            edgecolors=mpl.colors.to_rgba(layer.edge_color, layer.edge_alpha)
            if options.mesh
            else "face",
            linewidths=layer.line_width if options.mesh else 0.15,
            antialiaseds=options.mesh,
        )
        collection.set_array(values[patch.cell_ids] * layer.factor)
        axis.add_collection(collection)
        axis.set_xlim(-patch.width / 2, patch.width / 2)
        axis.set_ylim(-patch.width / 2, patch.width / 2)
        axis.set_aspect("equal")
        axis.set_axis_off()
        canvas.draw()
        return np.asarray(canvas.buffer_rgba()).copy()
    finally:
        art.clear()


def save_panel(pixels: np.ndarray, destination: Path, options: TissueOptions) -> None:
    """Export an exact square image without titles, colorbars, margins or projection."""
    art = Figure(figsize=(3, 3), dpi=options.dpi)
    FigureCanvasAgg(art)
    try:
        axis = art.add_axes((0, 0, 1, 1))
        axis.imshow(pixels, interpolation="none")
        axis.set_axis_off()
        for extension in options.formats:
            art.savefig(
                destination.with_suffix(f".{extension}"), dpi=options.dpi, pad_inches=0
            )
    finally:
        art.clear()


def render_tissue(
    data: PublicationData,
    source_study: str,
    source_run: str,
    ranges: dict[str, tuple[float, float]],
    directory: Path,
    stem: str,
    options: TissueOptions,
) -> dict[str, object]:
    """Write standalone map panels and one annotated three-field time grid."""
    source = data.source(source_study, source_run)
    samples = data.snapshots(source, options.times)
    patch = tile_lobules(samples.faces, samples.cell_ids, options.width)
    art = figure((2.5 * len(options.times) + 1.1, 8.5))
    grid = art.add_gridspec(
        4,
        len(options.times) + 1,
        height_ratios=[0.25, 1, 1, 1],
        width_ratios=[1] * len(options.times) + [0.07],
        hspace=0.03,
        wspace=0.03,
    )
    heading = art.add_subplot(grid[0, :-1])
    heading.set_axis_off()
    title = heading.text(
        0.5,
        0.9,
        f"B  Linear increase | D{options.dose} | {options.width:g} x {options.width:g} lobule widths",
        ha="center",
        va="bottom",
        fontsize=14,
        fontweight="bold",
        transform=heading.transAxes,
    )
    add_pattern_icon(title, 1)
    heading.annotate(
        "",
        xy=(0.97, 0.3),
        xytext=(0.1, 0.3),
        xycoords="axes fraction",
        arrowprops={"arrowstyle": "-|>", "lw": 1.4, "color": "black"},
    )
    heading.text(
        0,
        0.3,
        "Time",
        va="center",
        fontsize=12,
        fontweight="bold",
        transform=heading.transAxes,
    )
    try:
        for row, layer in enumerate(TISSUE_LAYERS):
            field = SUMMARY_FIELDS.index(layer.sid)
            limits = display_limits(layer, ranges)
            for column, time in enumerate(options.times):
                pixels = patch_image(
                    patch, samples.frame(time)[field], layer, limits, options
                )
                save_panel(pixels, directory / panel_stem(stem, layer, time), options)
                axis = art.add_subplot(grid[row + 1, column])
                half = patch.width / 2
                axis.imshow(
                    pixels, extent=(-half, half, -half, half), interpolation="none"
                )
                axis.set_axis_off()
                if row == 0:
                    axis.set_title(f"{time / 3600:g} h", fontsize=12, fontweight="bold")
                if column == 0:
                    axis.text(
                        -0.05,
                        0.5,
                        layer.title.replace(" [mM]", ""),
                        transform=axis.transAxes,
                        rotation=90,
                        ha="right",
                        va="center",
                        fontsize=12,
                        fontweight="bold",
                    )
            bar = art.add_subplot(grid[row + 1, -1])
            if layer.sid == "rr_necrosis":
                bar.set_axis_off()
                handles = [
                    Patch(facecolor="black", edgecolor="black", label="necrosis"),
                    Patch(facecolor="white", edgecolor="black", label="no necrosis"),
                ]
                if options.vessel_color == "lightgray":
                    handles.append(
                        Patch(
                            facecolor="lightgray",
                            edgecolor="black",
                            label="vessel lumen",
                        )
                    )
                art.legend(
                    handles=handles,
                    loc="center left",
                    bbox_to_anchor=(0, 0.5),
                    bbox_transform=bar.transAxes,
                    frameon=False,
                    fontsize=11,
                )
            else:
                shared_colorbar(
                    art,
                    bar,
                    ScalarMappable(norm=Normalize(*limits), cmap=layer.colormap),
                    layer,
                    scale=1,
                    vertical=True,
                )
        save_publication(
            art,
            directory / stem,
            RenderOptions(formats=options.formats, dpi=options.dpi),
        )
    finally:
        art.clear()
    return {
        "lobule_count": len(patch.centers),
        "lobule_centers_in_lobule_widths": patch.centers.tolist(),
        "source_lobule_width": patch.source_width,
        "source_coordinate_units": "unchanged; no independent length-unit metadata",
        "square_area_in_lobule_widths_squared": patch.nominal_area,
        "native_tissue_area_in_lobule_widths_squared": patch.tissue_area,
        "native_vessel_void_fraction": max(
            0.0, 1 - patch.tissue_area / patch.nominal_area
        ),
        "clipped_triangles": len(patch.faces),
    }


def panel_stem(stem: str, layer: Layer, time: float) -> str:
    """Keep field and exact time in unambiguous, portable panel filenames."""
    token = (
        repr(float(time))
        .removesuffix(".0")
        .replace(".", "p")
        .replace("-", "m")
        .replace("+", "")
    )
    return f"{stem}-{FIELD_NAMES[layer.sid]}-t{token}s"


def create_tissue_panels(
    data_root: Path, options: TissueOptions | None = None, *, force: bool = False
) -> Path:
    """Publish additional flat tissue images from the existing local SPT caches."""
    options = options or TissueOptions()
    options.validate()
    data = PublicationData(data_root)
    sources = simulation_cohort(data)
    selected = [
        source
        for source in sources
        if source.run.condition
        and source.run.condition.pattern_key == 1
        and source.dose == options.dose - 1
    ]
    if len(selected) != 1:
        raise ValueError(
            "Exactly one linear-increase simulation is required for the selected dose"
        )
    source = selected[0]
    if min(options.times) < source.start or max(options.times) > source.end:
        raise ValueError(
            "Tissue snapshot times are outside the native simulation history"
        )
    ranges = data.ranges(sources)
    paths = data.dependencies(sources, [(source.run.study, source.run.run_id)])
    stem = f"linear-increase-tissue-d{options.dose}"
    title = f"Linear increase | D{options.dose} | {options.width:g} x {options.width:g} lobule widths"
    titles = {stem: title}
    orders = {stem: 2000}
    for row, layer in enumerate(TISSUE_LAYERS):
        for column, time in enumerate(options.times):
            name = panel_stem(stem, layer, time)
            titles[name] = f"{title} | {layer.title} at {time / 3600:g} h"
            orders[name] = 2001 + row * len(options.times) + column
    configuration: dict[str, object] = {
        "version": TISSUE_VERSION,
        "title": title,
        "study": source.run.study,
        "runs": [source.run.run_id],
        "artifact_titles": titles,
        "artifact_order": orders,
        "options": json.loads(
            json.dumps(
                {
                    **asdict(options),
                    "layers": [asdict(layer) for layer in TISSUE_LAYERS],
                }
            )
        ),
        "geometry_view": {
            "projection": "orthographic XY",
            "tilt_from_top_degrees": 0,
            "thickness": "not displayed",
        },
        "width_definition": "native flat-to-flat lobule width",
        "protein_field": "rr_protein displayed as CYP2E1",
        "scalar_ranges": {
            layer.sid: list(display_limits(layer, ranges)) for layer in TISSUE_LAYERS
        },
        "composition": "identical native lobules translated on a hexagonal lattice; visualization only",
        "native_cell_values": "unchanged, including clipped border cells; no averaging",
        "vessel_openings": "native lumens retained as background; no fabricated scalar data",
    }
    directory = data.root / "publication" / stem
    with run_lock(directory):
        cached = current_manifest(directory, configuration, data.root)
        if cached and not force:
            generation = generation_path(directory, cached)
            print(f"Tissue panels D{options.dose}: cached", flush=True)
        else:
            records = [fingerprint(path) for path in paths]
            with tempfile.TemporaryDirectory(
                prefix=".partial-", dir=directory
            ) as temporary:
                staging = Path(temporary)
                with mpl.rc_context({"svg.fonttype": "none", "pdf.fonttype": 42}):
                    layout = render_tissue(
                        data,
                        source.run.study,
                        source.run.run_id,
                        ranges,
                        staging,
                        stem,
                        options,
                    )
                description = (
                    f"Linear-increase tissue panels, D{options.dose} ({source.run.run_id}). "
                    f"A square {options.width:g} native flat-to-flat lobule widths per side is completely covered by a hexagonal lattice. "
                    f"Native vessel lumens remain {options.vessel_color} and carry no scalar data. "
                    "CYP2E1 displays the model's rr_protein field; toxic compound T and necrosis use their original fields. "
                    "Continuous native cell values are interpolated in time; necrosis holds its preceding native state. "
                    "Copies retain identical cell values, including at clipped borders. This is a repeated-lobule visualization, not a new multi-lobule transport simulation. "
                    "The view is directly from above, without tilt, perspective or visible thickness. "
                    "Standalone panels are exact squares without labels or outer margins; the composite supplies times, shared scalar bars and a binary necrosis legend."
                )
                atomic_text(staging / "caption.md", description + "\n")
                atomic_json(
                    staging / "provenance.json",
                    {**configuration, "layout": layout, "raw_hdf_reads": 0},
                )
                generation = publish(directory, staging, configuration, paths, records)
            print(f"Tissue panels D{options.dose}: rendered", flush=True)
        for name in titles:
            export_figure(generation, data.root, name)
    return generation


def main(argv: list[str] | None = None) -> int:
    """Create publication-ready square tissue cutouts and refresh the HTML report."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, default=RESULTS_DIR / "spt")
    parser.add_argument("--dose", type=int, default=6)
    parser.add_argument(
        "--times-hours",
        type=float,
        nargs="+",
        default=[time / 3600 for time in SIMULATION_TIMES],
    )
    parser.add_argument("--width-lobules", type=float, default=2.5)
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument(
        "--formats",
        nargs="+",
        choices=("svg", "pdf", "png"),
        default=["svg", "pdf", "png"],
    )
    parser.add_argument(
        "--vessel-color", choices=("white", "lightgray"), default="white"
    )
    parser.add_argument("--mesh", action="store_true")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args(argv)
    try:
        options = TissueOptions(
            dose=args.dose,
            times=tuple(time * 3600 for time in args.times_hours),
            width=args.width_lobules,
            dpi=args.dpi,
            formats=tuple(args.formats),
            vessel_color=args.vessel_color,
            mesh=args.mesh,
        )
        generation = create_tissue_panels(args.data_root, options, force=args.force)
        atomic_json(
            args.data_root / "tissue-panel-report.json",
            {
                "status": "passed",
                "generation": str(generation),
                "options": asdict(options),
            },
        )
        print(build_report(args.data_root))
    except (OSError, ValueError, KeyError) as error:
        print(f"Tissue panels failed: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
