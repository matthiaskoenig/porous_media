"""Simulation conditions and safe extraction of the study's generation settings."""

import ast
import math
from dataclasses import dataclass
from pathlib import Path

boundary_flows = [
    -1.60747e-09,
    -7.58061e-06,
    -1.41387e-05,
    -2.51064e-05,
    -4.15705e-05,
    -5.37477e-05,
    -6.33526e-05,
    -7.12447e-05,
]

pattern_idx2name = {
    0: "constant",
    1: "linear increase",
    2: "sharp pericentral",
    3: "linear decrease",
    4: "sharp periportal",
    5: "random",
}


@dataclass(frozen=True)
class SimulationCondition:
    """Identify a simulation without deriving its meaning from file ordering."""

    sim_id: str
    pattern_key: int
    boundary_flow_key: int
    boundary_flow: float


def simulation_conditions(
    flows: list[float] | None = None, patterns: list[int] | None = None
) -> list[SimulationCondition]:
    """Create stable IDs in the same pattern/flow order as the FEBio generator."""
    flows = boundary_flows if flows is None else flows
    patterns = list(pattern_idx2name) if patterns is None else patterns
    if not flows or any(not math.isfinite(flow) for flow in flows):
        raise ValueError("Boundary flows must be nonempty and finite")
    if len(set(flows)) != len(flows) or len(set(patterns)) != len(patterns):
        raise ValueError("Patterns and boundary flows must be unique")
    if not patterns or any(pattern not in pattern_idx2name for pattern in patterns):
        raise ValueError("Unknown or empty zonation patterns")
    return [
        SimulationCondition(f"sim{index:03}", pattern, flow_index, flow)
        for index, (pattern, flow_index, flow) in enumerate(
            (
                (pattern, flow_index, flow)
                for pattern in patterns
                for flow_index, flow in enumerate(flows)
            ),
            start=1,
        )
    ]


def generation_script(root: Path) -> Path:
    """Return the study's archived generation script, which is never executed."""
    return root / "simulations_sixth/material/feb_manipulation_v2.py"


def conditions_from_source(root: Path) -> tuple[list[SimulationCondition], str]:
    """Read literal generator settings or label the repository defaults explicitly.

    Arbitrary Python expressions and duplicate assignments are rejected. Raw
    generation scripts are parsed as data, never imported or executed.
    """
    script = generation_script(root)
    if not script.is_file():
        return simulation_conditions(), "repository defaults (unverified)"
    assignments: dict[str, object] = {}
    tree = ast.parse(script.read_text(encoding="utf-8"))
    loops = [
        node
        for node in ast.walk(tree)
        if isinstance(node, ast.For)
        and isinstance(node.iter, ast.Name)
        and node.iter.id == "zonation_pattern"
    ]
    if len(loops) != 1 or not any(
        isinstance(node, ast.For)
        and isinstance(node.iter, ast.Name)
        and node.iter.id == "boundary_flow"
        for node in loops[0].body
    ):
        raise ValueError(
            "Unrecognized generator ordering; expected pattern then flow loops"
        )
    for node in ast.walk(tree):
        if not isinstance(node, ast.Assign):
            continue
        for target in node.targets:
            if isinstance(target, ast.Name) and target.id in {
                "boundary_flow",
                "zonation_pattern",
            }:
                if target.id in assignments:
                    raise ValueError(f"Ambiguous generator assignment: {target.id}")
                assignments[target.id] = ast.literal_eval(node.value)
    flows = assignments.get("boundary_flow")
    patterns = assignments.get("zonation_pattern")
    if not isinstance(flows, list) or not isinstance(patterns, list):
        raise ValueError("Generator must contain literal flow and pattern lists")
    if any(type(value) not in (int, float) for value in flows) or any(
        type(value) is not int for value in patterns
    ):
        raise ValueError("Generator settings must be numeric literals")
    return simulation_conditions(
        [float(value) for value in flows], [int(value) for value in patterns]
    ), str(script)
