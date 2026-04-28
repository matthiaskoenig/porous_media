# SPT analysis

These scripts generate the visualizations for the SPT convergence study and analysis.

## Run analysis
- run simulations on cluster
- copy results in data folder
- update data path in [`__init__.py`](__init__.py)
- update boundary flows and run [`information.py`](information.py)
- process data into xdmf, lobule reconstruction, calculate xarrays [`spt_data_processing.py`](data_processing.py)
- run [`plots_convergence.py`](plots_convergence.py)
- run [`plots_convergence_geometry.py`](plots_convergence_geometry.py)
- run [`plots_simulations.py`](plots_simulations.py)
- run [`plots_simulations_geometry.py`](plots_simulations_geometry.py)

## 2D Geometry visualization
- zonation pattern plot
- videos (all simulations)
- snapshots of selected timepoints (all simulations)
