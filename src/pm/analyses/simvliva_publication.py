"""
Visualization of Simliva results.
"""
from pathlib import Path
from typing import Dict, List

import numpy as np

from pm import DATA_DIR
from pm.console import console
from pm.mesh.mesh_tools import MeshTimepoint

vtk_dirs: Dict[str, Path] = {
    'sim_T277': DATA_DIR / "simliva" / "005_T_277_15K_P0__0Pa_t_24h" / "vtk",
    'sim_T310': DATA_DIR / "simliva" / "006_T_310_15K_P0__0Pa_t_24h" / "vtk",
}

vtks: Dict[str, List[Path]] = {k: sorted(list(p.glob("*.vtk"))) for k, p in vtk_dirs.items()}
# console.print(vtks)
for sim_key, vtk_paths in vtks.items():
    console.print(f"{sim_key}: {len(vtk_paths)}")



scalars_iri = {
    'necrosis': {"title": "Necrosis (0: alive, 1: death)", "cmap": "binary"},
    'ATP': {"title": "ATP (mM)", "cmap": "RdBu"},
    'GLC': {"title": "Glucose (mM)", "cmap": "RdBu"},
    'LAC': {"title": "Lactate (mM)", "cmap": "RdBu"},
    'O2': {"title": "Oxygen (mM)", "cmap": "RdBu"},
    'PYR': {"title": "Pyruvate (mM)", "cmap": "RdBu"},
    'ADP': {"title": "ADP (mM)", "cmap": "RdBu"},
    'NADH': {"title": "NADH (mM)", "cmap": "RdBu"},
    'NAD': {"title": "NAD (mM)", "cmap": "RdBu"},
    'ROS': {"title": "ROS (mM)", "cmap": "RdBu"},
    'ALT': {"title": "ALT (mM)", "cmap": "RdBu"},
    'AST': {"title": "AST (mM)", "cmap": "RdBu"},
}


def calculate_limits(vtks: Dict[str, List[Path]]):
    """Calculate the limits from given VTKs for scalars."""
    # Get limits of variables
    limits: Dict[str, List[float]] = {k: None for k in scalars_iri}

    for sim_key, vtk_paths in vtks.items():
        for vtk_path in vtk_paths:
            console.print(vtk_path)
            mesh_tp: MeshTimepoint = MeshTimepoint.from_vtk(vtk_path=vtk_path, show=False)
            for scalar in scalars_iri:
                data = mesh_tp.cell_data[scalar]

                # new min, max
                dmin = data.min()
                dmax = data.max()
                scalar_limits = limits[scalar]
                if scalar_limits is None:
                    scalar_limits = [dmin, dmax]
                else:
                    if scalar_limits[0] > dmin:
                        scalar_limits[0] = dmin
                    if scalar_limits[1] < dmax:
                        scalar_limits[1] = dmax

                limits[scalar] = scalar_limits
                # console.print(data)
                # console.print(scalar_limits)

    console.print(limits)
    return limits

limits = {
    'necrosis': [0.0, 1.0],
    'ATP': [0.1, 3.35109],
    'GLC': [0.0, 5.0],
    'LAC': [1.0, 4.164],
    'O2': [0.0, 1.2],
    'PYR': [0.680262, 2.77601],
    'ADP': [0.648908, 3.9],
    'NADH': [3.9, 4.0],
    'NAD': [0.0, 0.1],
    'ROS': [0.0, 0.826532],
    'ALT': [0.0, 1.0],
    'AST': [0.0, 0.7]
}
# limits = calculate_limits(vtks)

# merge limits
for scalar, lims in limits.items():
    scalars_iri[scalar]["clim"] = lims

console.print(scalars_iri)


from pm import BASE_DIR
output_path = BASE_DIR / "results" / "simliva_publication"


if __name__ == "__main__":

    # Create figures for all simulation timepoints of relevance (skip some results)
    for sim_key, vtk_paths in vtks.items():

        panels_path = output_path / "panels" / sim_key
        panels_path.mkdir(exist_ok=True, parents=True)

        # Create all figures for all variables
        from pm.visualization.pyvista_visualization import visualize_lobulus_vtk

        k = 0
        for vtk_path in vtk_paths:
            visualize_lobulus_vtk(
                vtk_path=vtk_path,
                scalars=scalars_iri,
                output_dir=panels_path,
            )
            k = k + 1
            if k>5:
                break


    # Create combined figures for variables and timepoints


