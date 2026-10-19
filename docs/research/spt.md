# SPT analysis

These historical notes describe the SPT convergence study and simulation
visualizations. The scripts require private research datasets, local path
configuration, and simulation results produced outside this package. They are
not a portable tutorial or part of the supported public API.

Use the [quickstart](../quickstart.md) for a reproducible example with bundled
data. See the [research archive](../research.md) for related historical notes.

## Analysis workflow

For researchers with access to the original simulation environment and datasets:

1. Run the simulations on the cluster and copy the results into the data folder.
2. Set `data_spt_dir` in
   [`__init__.py`](https://github.com/matthiaskoenig/porous_media/blob/main/src/porous_media/analyses/spt/__init__.py)
   to the local dataset directory.
3. Update the boundary flows and run
   [`information.py`](https://github.com/matthiaskoenig/porous_media/blob/main/src/porous_media/analyses/spt/information.py).
4. Review and run the relevant stages in
   [`data_processing.py`](https://github.com/matthiaskoenig/porous_media/blob/main/src/porous_media/analyses/spt/data_processing.py)
   to convert results to XDMF, reconstruct lobules, and calculate xarray datasets.
5. Generate the convergence plots with
   [`plots_convergence.py`](https://github.com/matthiaskoenig/porous_media/blob/main/src/porous_media/analyses/spt/plots_convergence.py)
   and
   [`plots_convergence_geometry.py`](https://github.com/matthiaskoenig/porous_media/blob/main/src/porous_media/analyses/spt/plots_convergence_geometry.py).
6. Generate the simulation plots with
   [`plots_simulations.py`](https://github.com/matthiaskoenig/porous_media/blob/main/src/porous_media/analyses/spt/plots_simulations.py)
   and
   [`plots_simulations_geometry.py`](https://github.com/matthiaskoenig/porous_media/blob/main/src/porous_media/analyses/spt/plots_simulations_geometry.py).

Review each script's dataset selections and output directories before running
it. The historical configuration reflects the original research environment.

## 2D geometry visualization

The study visualizations include:

- Zonation pattern plots.
- Videos for all simulations.
- Snapshots at selected time points for all simulations.
