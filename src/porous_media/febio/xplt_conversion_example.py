from pathlib import Path

import porous_media.febio.febpy.febio_vtk as fevtk
import porous_media.febio.febpy.xplt as xplt


# read the xplt file
test_data_path = Path(__file__).parent / "lobule_zonation_pattern_sim001.xplt"
reader = xplt.xplt(test_data_path)
reader.readAllStates()

# Convert
wall = fevtk.DomainToVTKConverter("wall", reader)
fluid = fevtk.DomainToVTKConverter("fluid", reader)

# create a vtk grid with data from time step 4
wall_vtk = wall.set_data_from_xplt(reader, 4)

_ = wall_vtk.set_active_scalars("position")
wall_vtk.plot()
# Continue working with the vtk object

# Save as a vtk file for external use
wall_vtk.save("test_4.vtk")
