using GridapGmsh
using GridapGmsh: gmsh, GmshDiscreteModel

function add_square!(gmsh, s, x0, y0)
    factory = gmsh.model.geo

    p1 = factory.addPoint(x0, y0, 0, lc)
    p2 = factory.addPoint(x0+s, y0, 0, lc)
    p3 = factory.addPoint(x0+s, y0+s, 0, lc)
    p4 = factory.addPoint(x0, y0+s, 0, lc)

    l1 = factory.addLine(p1, p2)
    l2 = factory.addLine(p2, p3)
    l3 = factory.addLine(p3, p4)
    l4 = factory.addLine(p4, p1)

    cl1 = factory.addCurveLoop([l1, l2, l3, l4])

    return cl1
end

p = 1.0
N = 10
lc = 1.0

gmsh.initialize()
gmsh.model.add("subcrit-assy")

fuelTags = Int32[]
sourceTags = Int32[]

fuel_pos = [(5,6)]
source_pos = [(8,3)]

for i in 1:N, j in 1:N
    if (i, j) in fuel_pos
        x0 = (i - 1) * p
        y0 = (j - 1) * p
        squareTag = add_square!(gmsh, p, x0, y0)
        push!(fuelTags, squareTag)
    elseif (i, j) in source_pos
        x0 = (i - 1) * p
        y0 = (j - 1) * p
        squareTag = add_square!(gmsh, p, x0, y0)
        push!(sourceTags, squareTag)
    end
end

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

h2oTag = factory.addPlaneSurface(vcat(cl1, fuelTags, sourceTags))
fuelTag = factory.addPlaneSurface(fuelTags)
sourceTag = factory.addPlaneSurface(sourceTags)

pg1 = gmsh.model.geo.addPhysicalGroup(2, [fuelTag])
pg2 = gmsh.model.geo.addPhysicalGroup(2, [sourceTag])
pg4 = gmsh.model.geo.addPhysicalGroup(2, [h2oTag])
pg5 = gmsh.model.geo.addPhysicalGroup(1, [l1])
pg6 = gmsh.model.geo.addPhysicalGroup(1, [l2])
pg7 = gmsh.model.geo.addPhysicalGroup(1, [l3])
pg8 = gmsh.model.geo.addPhysicalGroup(1, [l4])

gmsh.model.setPhysicalName(2, pg1, "fuel")
gmsh.model.setPhysicalName(2, pg2, "source")
gmsh.model.setPhysicalName(2, pg4, "water")
gmsh.model.setPhysicalName(1, pg5, "bottom")
gmsh.model.setPhysicalName(1, pg6, "right")
gmsh.model.setPhysicalName(1, pg7, "top")
gmsh.model.setPhysicalName(1, pg8, "left")


gmsh.model.geo.removeAllDuplicates()

gmsh.model.geo.synchronize()

gmsh.model.mesh.generate(2)

gmsh.write("subcrit-assy.msh")

if !("-nopopup" in ARGS)
    gmsh.fltk.run()
end

gmsh.finalize()

# move to json file format
using Gridap
mshfile = joinpath(@__DIR__,"subcrit-assy.msh")
model = GmshDiscreteModel(mshfile; renumber=true)
Gridap.Io.to_json_file(model, "subcrit-assy.json")