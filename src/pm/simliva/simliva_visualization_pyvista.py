"""Visualization with pyvista"""
from pm.console import console
import pyvista as pv

# global configuration
pv.global_theme.window_size = [1200, 1200]
pv.global_theme.background = 'white'
pv.global_theme.transparent_background = True
# pv.global_theme.cmap = 'RdBu'
pv.global_theme.colorbar_orientation = 'vertical'

pv.global_theme.font.family = 'arial'
pv.global_theme.font.size = 20
pv.global_theme.font.title_size = 40
pv.global_theme.font.label_size = 20
pv.global_theme.font.color = 'black'


# read the data
grid = pv.read('lobule_BCflux.t006.vtk')
dset = grid.cell_data
# dset.

console.print(grid)
console.print(grid.cell_data)

# deactivate active sets
# grid.set_active_tensors(None)
grid.set_active_scalars(None)
# grid.set_active_vectors(None)

# plot the data with an automatically created Plotter
# grid.plot(show_scalar_bar=False, show_axes=False)


p = pv.Plotter(
    window_size=(1200, 1200),
    title="TPM lobulus",
    # shape=(3, 2), border=False,
    off_screen=True
)

# TODO: object position (i.e. position in box)

# p.subplot(0, 0)
grid.set_active_scalars(name="rr_(S)")


actor = p.add_mesh(
    grid,
    show_edges=True,
    render_points_as_spheres=True,
    point_size=5,
    show_vertices=True,
    line_width=1.0,
    cmap='RdBu',
    show_scalar_bar=False,
    # specular=0.5, specular_power=15
)

p.add_scalar_bar(
    title="Substrate S [mM]",
    n_labels=5,
    bold=True,
    # height=0.6,
    width=0.6,
    vertical=False,
    position_y=0,
    position_x=0.2,
    mapper=actor.mapper,
)

# TODO: set global min and max values on colormap (clim on Dat
# p.update_scalar_bar_range(clim=[0, 1], name="Substrate S [mM]")

# create all figures for data


# TODO: camera position
# Camera position to zoom to face
p.camera_position = (0, 3E-4, 1E-3)
p.camera.zoom(1.0)
print(p.camera)

p.show(
    # cpos=top_view,
    screenshot='test_pyvista.png'
)

