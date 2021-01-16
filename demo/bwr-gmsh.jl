using GridapGmsh
using GridapGmsh: gmsh, GmshDiscreteModel

function add_pin!(gmsh, o, r, t, lc)
    factory = gmsh.model.geo

    # inner and outer circle
    l1 = add_circle!(gmsh, o, r, lc)
    l2 = add_circle!(gmsh, o, r + t, lc)

    s1 = factory.addPlaneSurface([l1])
    s2 = factory.addPlaneSurface([l2, l1]) # l1 is a hole

    return s1, s2
end

function add_circle!(gmsh, o, r, lc)
    factory = gmsh.model.geo
    ox, oy = o

    p1 = factory.addPoint(ox, oy, 0, lc)  # center - origin
    p2 = factory.addPoint(ox + r, oy, 0, lc)  # right
    p3 = factory.addPoint(ox, oy + r, 0, lc)  # up
    p4 = factory.addPoint(ox - r, oy, 0, lc)  # left
    p5 = factory.addPoint(ox, oy - r, 0, lc)  # down

    c1 = factory.addCircleArc(p2, p1, p3)
    c2 = factory.addCircleArc(p3, p1, p4)
    c3 = factory.addCircleArc(p4, p1, p5)
    c4 = factory.addCircleArc(p5, p1, p2)

    l1 = factory.addCurveLoop([c1, c2, c3, c4])

    return l1
end

function add_square!(gmsh, s)
    factory = gmsh.model.geo

    p1 = factory.addPoint(0, 0, 0, lc)
    p2 = factory.addPoint(s, 0, 0, lc)
    p3 = factory.addPoint(s, s, 0, lc)
    p4 = factory.addPoint(0, s, 0, lc)

    l1 = factory.addLine(p1, p2)
    l2 = factory.addLine(p2, p3)
    l3 = factory.addLine(p3, p4)
    l4 = factory.addLine(p4, p1)

    cl1 = factory.addCurveLoop([l1, l2, l3, l4])

    return cl1
end

N = 4

# in cm
p = 1.6        # pitch
ri = 0.5       # internal radius
t = 0.1        # wall thickness
ro = ri + t    # external radius
lc = 0.1

gmsh.initialize()
gmsh.model.add("bwr")

# tags
pinTags = Int32[]
cladTags = Int32[]
gdPinTags = Int32[]

# gd pins positions
GD_pos = [(2, 3), (3, 2)]

for i in 1:N, j in 1:N
    xo = p / 2 + (i - 1) * p
    yo = p / 2 + (j - 1) * p

    pinTag, cladTag = add_pin!(gmsh, (xo, yo), ri, t, lc)
    push!(cladTags, cladTag)
    if (i, j) in GD_pos
        push!(gdPinTags, pinTag)
    else
        push!(pinTags, pinTag)
    end
end

#! TODO: h2oTag = add_reflector!(gmsh, 4p)
s = N * p
factory = gmsh.model.geo

p1 = factory.addPoint(0, 0, 0, lc)
p2 = factory.addPoint(s, 0, 0, lc)
p3 = factory.addPoint(s, s, 0, lc)
p4 = factory.addPoint(0, s, 0, lc)

l1 = factory.addLine(p1, p2)
l2 = factory.addLine(p2, p3)
l3 = factory.addLine(p3, p4)
l4 = factory.addLine(p4, p1)

cl1 = factory.addCurveLoop([l1, l2, l3, l4])

h2oTag = factory.addPlaneSurface(vcat(cl1, pinTags, cladTags, gdPinTags))

pg1 = gmsh.model.geo.addPhysicalGroup(2, pinTags)
pg2 = gmsh.model.geo.addPhysicalGroup(2, cladTags)
pg3 = gmsh.model.geo.addPhysicalGroup(2, gdPinTags)
pg4 = gmsh.model.geo.addPhysicalGroup(2, [h2oTag])
pg5 = gmsh.model.geo.addPhysicalGroup(1, [l1])
pg6 = gmsh.model.geo.addPhysicalGroup(1, [l2])
pg7 = gmsh.model.geo.addPhysicalGroup(1, [l3])
pg8 = gmsh.model.geo.addPhysicalGroup(1, [l4])

gmsh.model.setPhysicalName(2, pg1, "pin")
gmsh.model.setPhysicalName(2, pg2, "cladding")
gmsh.model.setPhysicalName(2, pg3, "pin-gd")
gmsh.model.setPhysicalName(2, pg4, "water")
gmsh.model.setPhysicalName(1, pg5, "bottom")
gmsh.model.setPhysicalName(1, pg6, "right")
gmsh.model.setPhysicalName(1, pg7, "top")
gmsh.model.setPhysicalName(1, pg8, "left")

# removemos puntos de la geometria duplicados (origenes por ejemplo)
gmsh.model.geo.removeAllDuplicates()

gmsh.model.geo.synchronize()

gmsh.model.mesh.generate(2)

gmsh.write("bwr.msh")

if !("-nopopup" in ARGS)
    gmsh.fltk.run()
end

gmsh.finalize()

# move to json file format
using Gridap
mshfile = joinpath(@__DIR__,"bwr.msh")
model = GmshDiscreteModel(mshfile; renumber=true)
Gridap.Io.to_json_file(model, "bwr.json")
