# Information about libraries, tools and formats

# Working with meshes
- meshes are stored as a collection of `points` and `cells`, which are objects made from `points`
- data can be associated either with the `points` or the cells of the mesh

## Formats:
- **XDMF**: The XDMF format supports time series with a shared mesh. You can write times series data using meshio with
- FEBioStudio can convert some other formats to the FEBio input specification. For instance, `NIKE3D [54]` and `Abaqus` input files can be imported in FEBioStudio and can be exported as a FEBio input file.

## meshio
https://github.com/nschloe/meshio
There are various mesh formats available for representing unstructured meshes. meshio can read and write all of the following and smoothly converts between them:

    Abaqus (.inp), ANSYS msh (.msh), AVS-UCD (.avs), CGNS (.cgns), DOLFIN XML (.xml), Exodus (.e, .exo), FLAC3D (.f3grid), H5M (.h5m), Kratos/MDPA (.mdpa), Medit (.mesh, .meshb), MED/Salome (.med), Nastran (bulk data, .bdf, .fem, .nas), Netgen (.vol, .vol.gz), Neuroglancer precomputed format, Gmsh (format versions 2.2, 4.0, and 4.1, .msh), OBJ (.obj), OFF (.off), PERMAS (.post, .post.gz, .dato, .dato.gz), PLY (.ply), STL (.stl), Tecplot .dat, TetGen .node/.ele, SVG (2D output only) (.svg), SU2 (.su2), UGRID (.ugrid), VTK (.vtk), VTU (.vtu), WKT (TIN) (.wkt), XDMF (.xdmf, .xmf).

## Gmsh 
A three-dimensional finite element mesh generator with built-in pre- and post-processing facilities
https://gmsh.info/
- GUI and python interface; 

## trimesh
https://trimesh.org/
- GUI Tool

Trimesh is a pure Python 3.7+ library for loading and using triangular meshes with an emphasis on watertight surfaces. The goal of the library is to provide a full featured and well tested Trimesh object which allows for easy manipulation and analysis, in the style of the Polygon object in the Shapely library.
- MeshLab offers a series of automatic, semi-manual and interactive filters to remove those geometric element generally considered “wrong” by most software and algorithms. It is possible to removing topological errors, duplicated and unreferenced vertices, small components, degenerated or intersecting faces, and many more geometrical and topological singularities. Using different automatic and interactive selection methods, is then possible to isolate and remove unwanted areas of your meshes and point clouds.
- A common need when processing a 3D model is to reduce its geometric complexity, creating a geometry with the same shape but with less triangles (or points). MeshLab offers different ways to simplify (decimate) triangulated surfaces, able to preserve geometrical detail and texture mapping, or to selectively reduce the number of points in a pointcloud. In other cases, the user may want to increase the number of triangles (or points): MeshLab also provides different subdivision schemes, remeshing and resampling filters to increase geometric complexity of 3D models, or to optimize point distribution and triangulation quality.
## shapely

## meshlab
MeshLab the open source system for processing and editing 3D triangular meshes.
It provides a set of tools for editing, cleaning, healing, inspecting, rendering, texturing and converting meshes. It offers features for processing raw data produced by 3D digitization tools/devices and for preparing models for 3D printing.


## pymesh
https://pymesh.readthedocs.io/en/latest/
PyMesh is a rapid prototyping platform focused on geometry processing. It provides a set of common mesh processing functionalities and interfaces with a number of state-of-the-art open source packages to combine their power seamlessly under a single developing environment.

## meshplex





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

# Algorithms
## scipy.spatial
*Delauny* triangulation
- scipy.spatial: https://docs.scipy.org/doc/scipy/reference/spatial.html; 

- Delaunay triangulation, convex hulls, and Voronoi diagrams
- Spatial transformations: https://docs.scipy.org/doc/scipy/reference/spatial.transform.html#module-scipy.spatial.transform

- Morphometrics

# Finite element simulation
https://fenicsproject.org/
