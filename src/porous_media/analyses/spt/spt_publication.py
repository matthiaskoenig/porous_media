"""Visualization of SPT results."""
from pathlib import Path
from typing import Dict, List

import numpy as np

from porous_media import DATA_DIR
from porous_media.console import console
from porous_media.mesh.mesh_tools import MeshTimepoint
from porous_media.visualization.pyvista_visualization import visualize_lobulus_vtk


# FIXME: range(1, 26),
vtk_dirs: Dict[str, Path] = {
    f"sim__{k:03d}": DATA_DIR / "spt" / f"sim__{k:03d}" for k in range(1, 2)
}

changes = {
    "sim001": {"boundary_pressure": 236.17, "boundary_flow": -3.376671607824531e-06},
    "sim002": {"boundary_pressure": 236.17, "boundary_flow": -6.753343215649062e-06},
    "sim003": {"boundary_pressure": 236.17, "boundary_flow": -1.0130014823473593e-05},
    "sim004": {"boundary_pressure": 236.17, "boundary_flow": -1.3506686431298124e-05},
    "sim005": {"boundary_pressure": 236.17, "boundary_flow": -1.6883358039122653e-05},
    "sim006": {"boundary_pressure": 343.08, "boundary_flow": -6.753343215649062e-06},
    "sim007": {"boundary_pressure": 343.08, "boundary_flow": -1.3506686431298124e-05},
    "sim008": {"boundary_pressure": 343.08, "boundary_flow": -2.0260029646947185e-05},
    "sim009": {"boundary_pressure": 343.08, "boundary_flow": -2.7013372862596247e-05},
    "sim010": {"boundary_pressure": 343.08, "boundary_flow": -3.3766716078245306e-05},
    "sim011": {"boundary_pressure": 449.99, "boundary_flow": -1.0130014823473594e-05},
    "sim012": {"boundary_pressure": 449.99, "boundary_flow": -2.026002964694719e-05},
    "sim013": {"boundary_pressure": 449.99, "boundary_flow": -3.0390044470420787e-05},
    "sim014": {"boundary_pressure": 449.99, "boundary_flow": -4.052005929389438e-05},
    "sim015": {"boundary_pressure": 449.99, "boundary_flow": -5.065007411736797e-05},
    "sim016": {"boundary_pressure": 556.9, "boundary_flow": -1.3506686431298124e-05},
    "sim017": {"boundary_pressure": 556.9, "boundary_flow": -2.7013372862596247e-05},
    "sim018": {"boundary_pressure": 556.9, "boundary_flow": -4.052005929389437e-05},
    "sim019": {"boundary_pressure": 556.9, "boundary_flow": -5.4026745725192495e-05},
    "sim020": {"boundary_pressure": 556.9, "boundary_flow": -6.753343215649061e-05},
    "sim021": {"boundary_pressure": 663.81, "boundary_flow": -1.6883358039122653e-05},
    "sim022": {"boundary_pressure": 663.81, "boundary_flow": -3.3766716078245306e-05},
    "sim023": {"boundary_pressure": 663.81, "boundary_flow": -5.065007411736797e-05},
    "sim024": {"boundary_pressure": 663.81, "boundary_flow": -6.753343215649061e-05},
    "sim025": {"boundary_pressure": 663.81, "boundary_flow": -8.441679019561327e-05},
}

S_inflow = np.linspace(0.2, 1.0, num=5)
Qli_neu = np.linspace(0.2, 1.0, num=5)
#
# # all vtk
# console.rule(title="Number of VTKs", align="left", style="white")
# vtks_all: Dict[str, List[Path]] = {k: sorted(list(p.glob("*.vtk"))) for k, p in vtk_dirs.items()}
# for sim_key, vtk_paths in vtks_all.items():
#     console.print(f"{sim_key}: {len(vtk_paths)}")
#
# times_all: Dict[str, List[float]] = {}
# for sim_key, vtk_paths in vtks_all.items():
#     times = []
#     for p in vtk_paths:
#         with open(p, "r") as f_vtk:
#             f_vtk.readline()  # skip first line
#             line = f_vtk.readline()
#             t = float(line.strip().split(" ")[-1])
#             times.append(t)
#     times_all[sim_key] = times
#
# tends =
# for sim_key, times in times_all.items():
#     console.print(f"{sim_key}: time: [{times[0]}, {times[-1]}]")
#
#
#
# # filter vtk by times
# # interpolation
# times_wanted = np.linspace(0, 600*60, 11)  # [s] (21 points in 600 min) # static image 28800; 480 [min]
# # times_wanted = np.linspace(0, 600*60, 201)  # [s] (21 points in 600 min) # gifs
#
#
# vtks: Dict[str, List[Path]] = {}
# times: Dict[str, List[float]] = {}
# for sim_key, vtk_paths in vtks_all.items():
#     tps = times_all[sim_key]
#     paths_wanted = []
#     times_actual = []
#
#     for t_wanted in times_wanted:
#         for k, t in enumerate(tps):
#             if t >= t_wanted:
#                 paths_wanted.append(vtk_paths[k])
#                 times_actual.append(t)
#                 break
#
#     vtks[sim_key] = paths_wanted
#     times[sim_key] = times_actual
#
#
# # filter vtk by index
# console.rule(title="Number of filtered VTKs", align="left", style="white")
# console.print(f"{times_wanted=}")
# for sim_key, vtk_paths in vtks.items():
#     console.print(f"{sim_key}: {len(vtk_paths)}, times: {times[sim_key]}")
#
# scalars_iri = {
#     'necrosis': {"title": "Necrosis (0: alive, 1: death)", "cmap": "binary"},
#     'ATP': {"title": "ATP (mM)", "cmap": "RdBu"},
#     'GLC': {"title": "Glucose (mM)", "cmap": "RdBu"},
#     'LAC': {"title": "Lactate (mM)", "cmap": "RdBu"},
#     'O2': {"title": "Oxygen (mM)", "cmap": "RdBu"},
#     'PYR': {"title": "Pyruvate (mM)", "cmap": "RdBu"},
#     'ADP': {"title": "ADP (mM)", "cmap": "RdBu"},
#     'NADH': {"title": "NADH (mM)", "cmap": "RdBu"},
#     'NAD': {"title": "NAD (mM)", "cmap": "RdBu"},
#     'ROS': {"title": "ROS (mM)", "cmap": "RdBu"},
#     'ALT': {"title": "ALT (mM)", "cmap": "RdBu"},
#     'AST': {"title": "AST (mM)", "cmap": "RdBu"},
# }
#
# def calculate_limits(vtks: Dict[str, List[Path]]):
#     """Calculate the limits from given VTKs for scalars."""
#     # Get limits of variables
#     limits: Dict[str, List[float]] = {k: None for k in scalars_iri}
#
#     for sim_key, vtk_paths in vtks.items():
#         for vtk_path in vtk_paths:
#             console.print(vtk_path)
#             mesh_tp: MeshTimepoint = MeshTimepoint.from_vtk(vtk_path=vtk_path, show=False)
#             for scalar in scalars_iri:
#                 data = mesh_tp.cell_data[scalar]
#
#                 # new min, max
#                 dmin = data.min()
#                 dmax = data.max()
#                 scalar_limits = limits[scalar]
#                 if scalar_limits is None:
#                     scalar_limits = [dmin, dmax]
#                 else:
#                     if scalar_limits[0] > dmin:
#                         scalar_limits[0] = dmin
#                     if scalar_limits[1] < dmax:
#                         scalar_limits[1] = dmax
#
#                 limits[scalar] = scalar_limits
#                 # console.print(data)
#                 # console.print(scalar_limits)
#
#     console.print(limits)
#     return limits
#
# limits = {
#     'necrosis': [0.0, 1.0],
#     'ATP': [0.1, 3.35109],
#     'GLC': [0.0, 5.0],
#     'LAC': [1.0, 4.164],
#     'O2': [0.0, 1.2],
#     'PYR': [0.680262, 2.77601],
#     'ADP': [0.648908, 3.9],
#     'NADH': [3.9, 4.0],
#     'NAD': [0.0, 0.1],
#     'ROS': [0.0, 0.826532],
#     'ALT': [0.0, 1.0],
#     'AST': [0.0, 0.7]
# }
# # limits = calculate_limits(vtks)
#
# # merge limits
# for scalar, lims in limits.items():
#     scalars_iri[scalar]["clim"] = lims
#
# console.rule(title="Scalars", align="left", style="white")
# console.print(scalars_iri)
#
#
#
# def visualize_panels(vtks: Dict[str, List[Path]], output_path: Path, scalars):
#     # Create figures for all simulation timepoints of relevance (skip some results)
#     for sim_key, vtk_paths in vtks.items():
#
#         panels_path = output_path / sim_key / "panels"
#         panels_path.mkdir(exist_ok=True, parents=True)
#
#         # Create all figures for all variables
#         for k, vtk_path in enumerate(vtk_paths):
#             show_scalar_bar = True
#             # if k == 0 or k == (len(vtk_paths)-1):
#             #     show_scalar_bar = True
#             visualize_lobulus_vtk(
#                 vtk_path=vtk_path,
#                 scalars=scalars,
#                 output_dir=panels_path,
#                 window_size=(600, 600),
#                 scalar_bar=show_scalar_bar,
#             )
#
#
# if __name__ == "__main__":
#
#     from porous_media import BASE_DIR
#     output_path = BASE_DIR / "results" / "simliva_publication"
#     visualize_panels(vtks=vtks, scalars=scalars_iri, output_path=output_path)
#
#     # Create combined figures for variables and timepoints
#     from porous_media.visualization.image_manipulation import merge_images
#
#     scalars_plot = ["GLC", "O2", "LAC", "ATP", "ADP", "ROS", "necrosis", "ALT", "AST"]
#
#     # static images
#     # for sim_key, vtk_paths in vtks.items():
#     #     row_dir: Path = output_path / sim_key / "rows"
#     #     row_dir.mkdir(parents=True, exist_ok=True)
#     #
#     #     rows: List[Path] = []
#     #     for p in vtk_paths:
#     #         row: List[Path] = []
#     #         for scalar in scalars_plot:
#     #             img_path = output_path / sim_key / "panels" / scalar / f"{p.stem}.png"
#     #             row.append(img_path)
#     #
#     #         row_image: Path = output_path / sim_key / "rows" / f"{p.stem}.png"
#     #         merge_images(paths=row, direction="horizontal", output_path=row_image)
#     #         rows.append(row_image)
#     #
#     #     console.print(rows)
#     #     image: Path = output_path / sim_key / f"{sim_key}.png"
#     #     merge_images(paths=rows, direction="vertical", output_path=image)
#
#     # gifs
#     for k in range(len(times_wanted)):
#         rows: List[Path] = []
#         for sim_key, vtk_paths in vtks.items():
#             row_dir: Path = output_path / sim_key / "rows"
#             row_dir.mkdir(parents=True, exist_ok=True)
#
#             vtk_path = vtk_paths[k]
#             row: List[Path] = []
#             for scalar in scalars_plot:
#                 img_path = output_path / sim_key / "panels" / scalar / f"{p.stem}.png"
#                 row.append(img_path)
#
#             row_image: Path = output_path / sim_key / "rows" / f"{p.stem}.png"
#             merge_images(paths=row, direction="horizontal", output_path=row_image)
#             rows.append(row_image)
#
#         console.print(rows)
#         image: Path = output_path / sim_key / "gifs" / f"{k:<03}.png"
#         merge_images(paths=rows, direction="vertical", output_path=image)
