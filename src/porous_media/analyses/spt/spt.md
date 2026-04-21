# SPT analysis

These scripts generate the visualizations for the SPT analysis and publication.

## Run analysis
- run simulations on cluster
- process results into xdmf using the `porous_media/src/analysis/spt_data_processing.py`
- copy data in the `porous_media/data/spt/<results_date>` folder
- update `results_date` in `porous_media/src/analysis/__init__.py`
- update and run `porous_media/src/analysis/spt_information.py` (update boundary flows)
- update and run `porous_media/src/analysis/spt_plots_analysis.py`
- update and run `porous_media/src/analysis/spt_plots_geometry.py`


## 2D Geometry visualization
- zonation pattern plot
- videos (all simulations)
- snapshots of selected timepoints (all simulations)
