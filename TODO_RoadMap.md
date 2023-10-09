# Analysis

## XPLT -> VTK conversion

## Geometry
- [ ] FIXME: calculate volume from geometry and add to the elements (compare to given volume); calculated cell data
- [ ] FIXME: calculate the inflow & outflow surfaces; use the cell type to get the inflow and outflow areas
- [ ] FIXME: calculate perfusion with vector data and geometry; visualization of perfusion

---------
## Vector fields and streamlines
- [ ] FIXME: visualize vector data: fluid-flux: vector field & streamlines
- [ ] FIXME: better combination of panels into figure; i.e. title and colorbar only on selected panels (i.e. different variants of the panels)

- [ ] FIXME: add calculated point and cell data to simulation; I.e. definition of formulas
      such as unit conversions [Pa] -> [mmHg];
      This should happen as early as possible, probably already when generating the XDMF. Additional data
      could also be injected later on.
      This is similar to the added geometry parameters for the mesh, reuse funcationality

---------
# Mesh generation
- [ ] Histology image with staining; mesh generation from image (spatialPCA) -> species comparison


