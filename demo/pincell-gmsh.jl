using GridapGmsh
using GridapGmsh: gmsh, GmshDiscreteModel

gmsh.initialize()
gmsh.model.add("pincell")

factory = gmsh.model.geo

# Dimensions in cm.
p = 1.6        # pitch
p_2 = p / 2
ri = 0.5       # internal radius
t = 0.1        # wall thickness
ro = ri + t    # external radius

lc = 4e-2

# Fuel circle.
factory.addPoint(p_2, p_2, 0, lc, 1)       # center
factory.addPoint(p_2 + ri, p_2, 0, lc, 2)  # right
factory.addPoint(p_2, p_2 + ri, 0, lc, 3)  # up
factory.addPoint(p_2 - ri, p_2, 0, lc, 4)  # left
factory.addPoint(p_2, p_2 - ri, 0, lc, 5)  # down
factory.addCircleArc(2, 1, 3, 1)
factory.addCircleArc(3, 1, 4, 2)
factory.addCircleArc(4, 1, 5, 3)
factory.addCircleArc(5, 1, 2, 4)
factory.addCurveLoop([1, 2, 3, 4], 1)

# Cladding outer circle.
# factory.addPoint(p_2, p_2, 0, lc, 1)     # Reuse the fuel center point.
factory.addPoint(p_2 + ro, p_2, 0, lc, 6)  # right
factory.addPoint(p_2, p_2 + ro, 0, lc, 7)  # up
factory.addPoint(p_2 - ro, p_2, 0, lc, 8)  # left
factory.addPoint(p_2, p_2 - ro, 0, lc, 9)  # down
factory.addCircleArc(6, 1, 7, 5)
factory.addCircleArc(7, 1, 8, 6)
factory.addCircleArc(8, 1, 9, 7)
factory.addCircleArc(9, 1, 6, 8)
factory.addCurveLoop([5, 6, 7, 8], 2)

# Square cell boundary.
factory.addPoint(0, 0, 0, lc, 10)
factory.addPoint(p, 0, 0, lc, 11)
factory.addPoint(p, p, 0, lc, 12)
factory.addPoint(0, p, 0, lc, 13)
factory.addLine(10, 11,  9)
factory.addLine(11, 12, 10)
factory.addLine(12, 13, 11)
factory.addLine(13, 10, 12)
factory.addCurveLoop([9, 10, 11, 12], 3)

factory.addPlaneSurface([1], 1)
# Region between the outer and inner circles; loop 1 is the fuel hole.
factory.addPlaneSurface([2, 1], 2)
# Moderator square; loops 1 and 2 are holes for the pin regions.
factory.addPlaneSurface([3, 2, 1], 3)

# Material surface groups.
factory.addPhysicalGroup(2, [1], 1)
factory.addPhysicalGroup(2, [2], 2)
factory.addPhysicalGroup(2, [3], 3)
GridapGmsh.gmsh.model.setPhysicalName(2, 1, "pin")
GridapGmsh.gmsh.model.setPhysicalName(2, 2, "cladding")
GridapGmsh.gmsh.model.setPhysicalName(2, 3, "water")

# Boundary curve groups.
factory.addPhysicalGroup(1,  [9], 4)
factory.addPhysicalGroup(1, [10], 5)
factory.addPhysicalGroup(1, [11], 6)
factory.addPhysicalGroup(1, [12], 7)
GridapGmsh.gmsh.model.setPhysicalName(1, 4, "bottom")
GridapGmsh.gmsh.model.setPhysicalName(1, 5, "right")
GridapGmsh.gmsh.model.setPhysicalName(1, 6, "top")
GridapGmsh.gmsh.model.setPhysicalName(1, 7, "left")

factory.synchronize()

gmsh.model.mesh.generate(2)

mshfile = joinpath(@__DIR__, "pincell.msh")
gmsh.write(mshfile)

if !("-nopopup" in ARGS)
    gmsh.fltk.run()
end

gmsh.finalize()

# Convert the mesh to Gridap's JSON format.
using Gridap
model = GmshDiscreteModel(mshfile; renumber=true)
Gridap.Io.to_json_file(model, joinpath(@__DIR__, "pincell.json"))
