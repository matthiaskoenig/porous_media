import gmsh
import numpy as np

gmsh.initialize()

"""
Lobule Data from JupyterNotebook by LM
"""
# diameter = Außenkreis, Kantenlänge = außenkreisradius; spitze zu spitze
diameter_lobule = 1  # [mm] ~ 1mm

#diameter_pv = 35  # [μm] ~ 35+-25μm, https://doi.org/10.1002/hep.510280206
#diameter_cv = 75  # [μm] 25-150 μm, mean 75μm; https://anatomypubs.onlinelibrary.wiley.com/doi/10.1002/ar.24560

z_height = 0.2  # [mm] ~ 2mm; https://anatomypubs.onlinelibrary.wiley.com/doi/10.1002/ar.24560

mesh_refinement = 0.0002 # [??] #5
gmsh.option.setNumber("Mesh.CharacteristicLengthMin", mesh_refinement);
gmsh.option.setNumber("Mesh.CharacteristicLengthMax", mesh_refinement);


#### Automatic Unit Conversion ###
diameter_lobule = diameter_lobule/1000 # from mm to m
z_height = z_height/1000 # from mm to m

print("diameter_lobule: ", diameter_lobule, " m")
print("z_height: ", z_height, " m")

# We first create one cube:
gmsh.model.occ.addRectangle(0, 0, 0, diameter_lobule/2*0.87, diameter_lobule/2, 1)


gmsh.model.geo.mesh.setTransfiniteSurface(1)
gmsh.model.occ.extrude([(2, 1)], 0.0, 0.0, z_height,[1], [1], recombine=True)

#gmsh.model.mesh.refine()
print(gmsh.model.occ.getEntities())
gmsh.model.occ.synchronize()


#gmsh.option.setNumber("Mesh.SubdivisionAlgorithm", 2)


gmsh.model.mesh.generate(3)

gmsh.option.setNumber("Mesh.RecombineAll", 3)
gmsh.option.setNumber("Mesh.MshFileVersion", 2.0)
#transfinite = True
#gmsh.option.setNumber("Mesh.SubdivisionAlgorithm", 2)

gmsh.write("box_occ_0002.msh")
gmsh.write("box_occ_0002.vtk")
