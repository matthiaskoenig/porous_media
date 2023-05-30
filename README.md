# Working with meshes
- meshes are stored as a collection of `points` and `cells`, which are objects made from `points`
- data can be associated either with the `points` or the cells of the mesh

## meshio
https://github.com/nschloe/meshio
There are various mesh formats available for representing unstructured meshes. meshio can read and write all of the following and smoothly converts between them

## Gmsh 
A three-dimensional finite element mesh generator with built-in pre- and post-processing facilities
https://gmsh.info/

## pymesh
https://pymesh.readthedocs.io/en/latest/
PyMesh is a rapid prototyping platform focused on geometry processing. It provides a set of common mesh processing functionalities and interfaces with a number of state-of-the-art open source packages to combine their power seamlessly under a single developing environment.

## pymeshlab

## mesh manipulation

# Visualization

## VTK
https://examples.vtk.org/site/
https://pypi.org/project/vtk/
https://vtk.org/

- python interface/library
- very low-level, complicated to write

## Visit 
https://visit-dav.github.io/visit-website/
Python scripting:
https://visit-sphinx-github-user-manual.readthedocs.io/en/develop/python_scripting/index.html

## pyvista
https://github.com/pyvista/pyvista

PyVista is:

- Pythonic VTK: a high-level API to the Visualization Toolkit (VTK)
- mesh data structures and filtering methods for spatial datasets
- 3D plotting made simple and built for large/complex data geometries


## mayavi
https://pypi.org/project/mayavi/
http://docs.enthought.com/mayavi/mayavi/

## vedo
https://vedo.embl.es/
Integration with many external libraries like FEniCS for PDE and finite element solutions


## vtkplotter
Interactive tools in notebooks
https://www.youtube.com/watch?v=raiIft8VeRU

## paraview

## 3D slicer

# Algorithms
## scipy.spatial
*Delauny* triangulation
- scipy.spatial: https://docs.scipy.org/doc/scipy/reference/spatial.html; 

- Delaunay triangulation, convex hulls, and Voronoi diagrams
- Spatial transformations: https://docs.scipy.org/doc/scipy/reference/spatial.transform.html#module-scipy.spatial.transform

- Morphometrics

# Finite element simulation
https://fenicsproject.org/
