import gmsh
import numpy as np

gmsh.initialize()

"""
Lobule Data from JupyterNotebook by LM
"""
# diameter = Außenkreis, Kantenlänge = außenkreisradius; spitze zu spitze
diameter_lobule = 1  # [mm] ~ 1mm

diameter_pv = 35  # [μm] ~ 35+-25μm, https://doi.org/10.1002/hep.510280206
diameter_cv = 75  # [μm] 25-150 μm, mean 75μm; https://anatomypubs.onlinelibrary.wiley.com/doi/10.1002/ar.24560

z_height = 0.2  # [mm] ~ 2mm; https://anatomypubs.onlinelibrary.wiley.com/doi/10.1002/ar.24560

mesh_refinement = 0.0000725 # [??] #5
gmsh.option.setNumber("Mesh.CharacteristicLengthMin", mesh_refinement);
gmsh.option.setNumber("Mesh.CharacteristicLengthMax", mesh_refinement);

#### Automatic Unit Conversion ###
diameter_lobule = diameter_lobule/1000 # from mm to m
z_height = z_height/1000 # from mm to m
diameter_pv = diameter_pv*10**(-6) # from μm to m
diameter_cv = diameter_cv*10**(-6) # from μm to m

# Calculate Corners of Hexagon
y_height = np.round(diameter_lobule/2, 8)
y_side = np.round(diameter_lobule/4, 8)
x_side = np.round(np.sqrt(3)*diameter_lobule/4, 8)

z_height = np.round(z_height, 8)
radius_pv = np.round(diameter_pv/2, 8)
radius_cv = np.round(diameter_cv/2, 8)

print("diameter_lobule: ", diameter_lobule, " m")
print("z_height: ", z_height, " m")
print("radius_pv: ", radius_pv, " m")
print("radius_cv: ", radius_cv, " m")

#Hack around for cut loop, since getValue is not working
Point_top_coord = [0, y_height, 0]
Point_bottom_coord = [0, -y_height, 0]
Point_left_top_coord = [-x_side, y_side, 0]
Point_left_bottom_coord = [-x_side, -y_side, 0]
Point_right_top_coord = [x_side, y_side, 0]
Point_right_bottom_coord = [x_side, -y_side, 0]


#gmsh.model.geo.addPoint(0, y_height, 0, point_top_coord)
point_top_coord = gmsh.model.occ.addPoint(0, y_height, 0)
point_bottom_coord = gmsh.model.occ.addPoint(0, -y_height, 0)
point_left_top_coord = gmsh.model.occ.addPoint(-x_side, y_side, 0)
point_left_bottom_coord = gmsh.model.occ.addPoint(-x_side, -y_side, 0)
point_right_top_coord = gmsh.model.occ.addPoint(x_side, y_side, 0)
point_right_bottom_coord = gmsh.model.occ.addPoint(x_side, -y_side, 0)

line1 = gmsh.model.occ.addLine(point_top_coord, point_left_top_coord)
line2 = gmsh.model.occ.addLine(point_left_top_coord, point_left_bottom_coord)
line3 = gmsh.model.occ.addLine(point_left_bottom_coord, point_bottom_coord)
line4 = gmsh.model.occ.addLine(point_bottom_coord, point_right_bottom_coord)
line5 = gmsh.model.occ.addLine(point_right_bottom_coord, point_right_top_coord)
line6 = gmsh.model.occ.addLine(point_right_top_coord, point_top_coord)

 # Define the line loop
line_loop = gmsh.model.occ.addCurveLoop([line1, line2, line3, line4, line5, line6])

# Define the plane surface
plane_surface = gmsh.model.occ.addPlaneSurface([line_loop])

# Define CV_circ
circle_cv = gmsh.model.occ.addDisk(0, 0, 0, radius_cv, radius_cv, 20)

# Define BooleanDifference for PV_circ
gmsh.model.occ.cut([(2, plane_surface)], [(2, circle_cv)])


# Subtract circles from the plane surface
for point in [Point_top_coord, Point_left_top_coord,
                  Point_left_bottom_coord, Point_bottom_coord,
                  Point_right_bottom_coord, Point_right_top_coord]:
    #ask in forum why this is not working
    #xyz = gmsh.model.getValue(0, point, [])
    #gmsh.model.geo.addPoint(xyz[0], xyz[1], 0.12, lc, 103)
    circle = gmsh.model.occ.addDisk(point[0], point[1], point[2], radius_pv, radius_pv)
    finale = gmsh.model.occ.cut([(2, plane_surface)], [(2, circle)])


gmsh.model.geo.mesh.setTransfiniteSurface(1)
gmsh.model.occ.extrude([(2, 1)], 0.0, 0.0, z_height,[1], [1], recombine=True)  #[4, 4], [0.5, 1]

#gmsh.model.mesh.refine()
print(gmsh.model.occ.getEntities())
gmsh.model.occ.synchronize()


#gmsh.option.setNumber("Mesh.SubdivisionAlgorithm", 2)


gmsh.model.mesh.generate(3)

gmsh.option.setNumber("Mesh.RecombineAll", 3)
gmsh.option.setNumber("Mesh.MshFileVersion", 2.0)
#transfinite = True
#gmsh.option.setNumber("Mesh.SubdivisionAlgorithm", 2)

gmsh.write("lobule_occ_0000725.msh")
gmsh.write("lobule_occ_0000725.vtk")
