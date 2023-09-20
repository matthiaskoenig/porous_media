# Analysis

- [ ] FIXME: visualize vector data: fluid-flux: vector field & streamlines

- [ ] FIXME: better combination of panels into figure; i.e. title and colorbar only on selected panels (i.e. different variants of the panels)

- [ ] FIXME: add calculated point and cell data to simulation; I.e. definition of formulas
      such as unit conversions [Pa] -> [mmHg];
      This should happen as early as possible, probably already when generating the XDMF. Additional data
      could also be injected later on.
- [ ] FIXME: calculate volume from geometry and add to the elements (compare to given volume); calculated cell data

=> functions for calculating on mesh and create DataFrame
- [ ] FIXME: visualization of lineplots:
  - [ ] necrosis_fraction ~ time
  - [ ] necrosis_fraction ~ position (video); histogram


- [ ] FIXME: calculate perfusion with vector data and geometry
- [ ] Create images/code for spatial analysis (position ~ values) -> create zonated mesh;

- [ ] FIXME: extend limits to full information for variables (how is this combined with the predefined information?);
      Full set of information required for variables;

# Mesh generation
- [ ] Histology image with staining; mesh generation from image
  - [ ] Create numpy array with data from image 
