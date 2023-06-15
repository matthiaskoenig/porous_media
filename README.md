# Working with meshes
- meshes are stored as a collection of `points` and `cells`, which are objects made from `points`
- data can be associated either with the `points` or the cells of the mesh

## meshio
https://github.com/nschloe/meshio
There are various mesh formats available for representing unstructured meshes. meshio can read and write all of the following and smoothly converts between them

There are various mesh formats available for representing unstructured meshes. meshio can read and write all of the following and smoothly converts between them:

    Abaqus (.inp), ANSYS msh (.msh), AVS-UCD (.avs), CGNS (.cgns), DOLFIN XML (.xml), Exodus (.e, .exo), FLAC3D (.f3grid), H5M (.h5m), Kratos/MDPA (.mdpa), Medit (.mesh, .meshb), MED/Salome (.med), Nastran (bulk data, .bdf, .fem, .nas), Netgen (.vol, .vol.gz), Neuroglancer precomputed format, Gmsh (format versions 2.2, 4.0, and 4.1, .msh), OBJ (.obj), OFF (.off), PERMAS (.post, .post.gz, .dato, .dato.gz), PLY (.ply), STL (.stl), Tecplot .dat, TetGen .node/.ele, SVG (2D output only) (.svg), SU2 (.su2), UGRID (.ugrid), VTK (.vtk), VTU (.vtu), WKT (TIN) (.wkt), XDMF (.xdmf, .xmf).


## Gmsh 
A three-dimensional finite element mesh generator with built-in pre- and post-processing facilities
https://gmsh.info/

## pymesh
https://pymesh.readthedocs.io/en/latest/
PyMesh is a rapid prototyping platform focused on geometry processing. It provides a set of common mesh processing functionalities and interfaces with a number of state-of-the-art open source packages to combine their power seamlessly under a single developing environment.


## FeBio mesh formats
FEBioStudio can convert some other formats to the FEBio input specification. For instance, `NIKE3D [54]` and `Abaqus` input files can be imported in FEBioStudio and can be exported as a FEBio input file.


## mesh manipulation

# Visualization
https://examples.vtk.org/site/

 I suggest using a dedicated vtk library, either bare vtk, the higher-level mayavi.mlab or my recently acquired favourite, pyvista. I'll elaborate on all this in detail.
 ## https://pypi.org/project/vtk/
https://vtk.org/

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


## paraview

## 3D slicer


# Finite element simulation
https://fenicsproject.org/
