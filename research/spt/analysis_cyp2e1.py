"""Export the CYP2E1 necrosis response over position and against its protein."""

from __future__ import annotations

import pandas as pd
import xarray as xr

from research.spt.conditions import pattern_name
from research.spt.plot_data import PlotSource

CYP2E1_PATTERN = "cyp2e1"
NECROSIS_FIELD = "rr_necrosis"
PROTEIN_FIELD = "rr_protein"
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
    if not selected:
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
