from pathlib import Path

from mayavi import mlab
from mayavi.modules.surface import Surface

file_name = str(Path(__file__).parent.parent.parent / "data" / "sim001_perf1_BC1.t025.vtk")

# create a new figure, grab the engine that's created with it
fig = mlab.figure()
engine = mlab.get_engine()

# open the vtk file, let mayavi figure it all out
vtk_file_reader = engine.open(file_name)

# plot surface corresponding to the data
surface = Surface()
engine.add_filter(surface, vtk_file_reader)

# block until figure is closed
mlab.show()
