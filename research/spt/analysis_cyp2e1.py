"""Export the CYP2E1 necrosis response over position and against its protein."""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import pandas as pd
import xarray as xr

from research.spt.conditions import pattern_name
from research.spt.fields import LABELS, display_factor
from research.spt.plot_data import PlotSource
from research.spt.plot_render import RenderOptions, figure, save_figure
from research.spt.plot_style import simulation_color

CYP2E1_PATTERN = "cyp2e1"
NECROSIS_FIELD = "rr_necrosis"
PROTEIN_FIELD = "rr_protein"
POSITION_LABEL = "Lobular position (0 periportal, 1 perivenous)"
IDENTITY_COLUMNS = (
    "run",
    "pattern",
    "pattern_name",
    "dose",
    "boundary_flow",
    "time",
    "position",
    "bin_left",
    "bin_right",
)


def _cyp2e1_sources(sources: list[PlotSource]) -> list[PlotSource]:
    """Select the measured-profile cohort by name, never by run identifier."""
    return [
        source
        for source in sources
        if source.run.condition is not None
        and pattern_name(source.run.condition.pattern_key) == CYP2E1_PATTERN
    ]


def _field_frame(profiles: xr.Dataset, field: str, prefix: str) -> pd.DataFrame:
    """Flatten one field of a single-time profile selection into tidy rows."""
    frame = profiles.sel(field=field).to_dataframe().reset_index()
    return frame.rename(
        columns={
            name: f"{prefix}_{name}"
            for name in ("mean", "std", "minimum", "maximum")
            if name in frame
        }
    )


def cyp2e1_frames(
    profiles: xr.Dataset, sources: list[PlotSource]
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return necrosis over position and necrosis paired with protein, at the last time.

    Both quantities are read at the same, final spatial time. Per-cell protein does
    not change, but the mesh deforms, so cells cross position bins and the binned
    protein profile of the first time is not the one of the last.
    """
    selected = _cyp2e1_sources(sources)
    stored = {str(field) for field in profiles.field.values}
    if not selected or not {NECROSIS_FIELD, PROTEIN_FIELD} <= stored:
        return pd.DataFrame(), pd.DataFrame()
    conditions = {
        source.run.run_id: source.run.condition
        for source in selected
        if source.run.condition is not None
    }
    time = float(profiles.time.values[-1])
    at_time = profiles.sel(run=list(conditions), time=time)
    necrosis = _field_frame(at_time, NECROSIS_FIELD, "necrosis")
    protein = _field_frame(at_time, PROTEIN_FIELD, "protein")
    necrosis["pattern"] = necrosis.run.map(
        {run: condition.pattern_key for run, condition in conditions.items()}
    )
    necrosis["pattern_name"] = necrosis.pattern.map(pattern_name)
    necrosis["dose"] = necrosis.run.map(
        {run: condition.boundary_flow_key + 1 for run, condition in conditions.items()}
    )
    necrosis["boundary_flow"] = necrosis.run.map(
        {run: condition.boundary_flow for run, condition in conditions.items()}
    )
    paired = necrosis.merge(
        protein[["run", "position", "protein_mean", "protein_std"]],
        on=["run", "position"],
        validate="one_to_one",
    )
    return (
        necrosis[
            [
                *IDENTITY_COLUMNS,
                "necrosis_mean",
                "necrosis_std",
                "necrosis_minimum",
                "necrosis_maximum",
                "cell_count",
                "volume",
            ]
        ],
        paired[
            [
                *IDENTITY_COLUMNS,
                "protein_mean",
                "protein_std",
                "necrosis_mean",
                "necrosis_std",
            ]
        ],
    )


def render_cyp2e1_figures(
    necrosis: pd.DataFrame,
    paired: pd.DataFrame,
    directory: Path,
    options: RenderOptions,
) -> list[Path]:
    """Draw necrosis over position and necrosis against protein, one line per flux."""
    if necrosis.empty or paired.empty:
        return []
    written: list[Path] = []
    hours = float(necrosis["time"].iloc[0]) / 3600
    scale = display_factor(NECROSIS_FIELD)
    panels = (
        ("cyp2e1-necrosis-position", necrosis, "position", POSITION_LABEL, 1.0),
        (
            "cyp2e1-necrosis-protein",
            paired,
            "protein_mean",
            f"CYP2E1 {LABELS[PROTEIN_FIELD]}",
            display_factor(PROTEIN_FIELD),
        ),
    )
    # Vector text keeps the flux labels searchable in the published SVG, exactly as
    # the other analysis figures are written.
    with mpl.rc_context({"svg.fonttype": "none", "pdf.fonttype": 42}):
        for name, frame, x_column, x_label, x_scale in panels:
            art = figure((6.4, 3.6))
            try:
                axis = art.add_subplot()
                for dose in sorted(frame["dose"].astype(int).unique().tolist()):
                    ordered = frame[frame["dose"] == dose].sort_values("position")
                    pattern = int(ordered["pattern"].astype(int).iloc[0])
                    axis.plot(
                        ordered[x_column] * x_scale,
                        ordered["necrosis_mean"] * scale,
                        color=simulation_color(pattern, dose),
                        marker="o",
                        markersize=3,
                        label=f"D{dose}",
                    )
                axis.set_xlabel(x_label)
                axis.set_ylabel(LABELS[NECROSIS_FIELD])
                axis.set_ylim(-0.02 * scale, 1.02 * scale)
                # The curves reach the top of the axes at every flux that causes
                # necrosis, so an inset legend would cover them.
                axis.legend(
                    title="Substrate flux",
                    fontsize="small",
                    title_fontsize="small",
                    loc="center left",
                    bbox_to_anchor=(1.02, 0.5),
                    borderaxespad=0.0,
                )
                art.suptitle(f"CYP2E1 necrosis at {hours:g} h")
                save_figure(art, directory / f"{name}.svg", options)
            finally:
                art.clear()
            written.extend(
                directory / f"{name}.{extension}" for extension in options.formats
            )
    return written
