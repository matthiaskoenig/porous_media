# SPT analysis

These scripts generate the visualizations for the SPT analysis and publication.

## Run analysis
- run simulations on cluster
- copy results in data folder
- update data paths in [`__init__.py`](__init__.py)
- update boundary flows and run [`spt_information.py`](spt_information.py)
- process data into xdmf, lobule reconstruction, calculate xarrays [`spt_data_processing.py`](./spt_data_processing.py)
- update and run [`spt_plots_analysis.py`](spt_plots_analysis.py)
- update and run `porous_media/src/analysis/spt_plots_geometry.py`


## 2D Geometry visualization
- zonation pattern plot
- videos (all simulations)
- snapshots of selected timepoints (all simulations)


# TODO
- rewrite scripts for working with external harddrive (decouple from package)
- parallelization
