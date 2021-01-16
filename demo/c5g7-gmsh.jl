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

N = 17

# in cm
p = 1.26
r = 0.54
lc = 0.3

gmsh.initialize()
gmsh.model.add("c5g7")

# tags
guidetube = Int32[]
UO₂ = Int32[]
MOX_43 = Int32[]
MOX_7 = Int32[]
MOX_87 = Int32[]
fissionchamber = Int32[]

GT_pos = [( 3, 6), ( 3, 9), (3, 12),
          ( 4, 4), ( 4, 14),
          ( 6, 3), ( 6, 6), ( 6, 9), ( 6, 12), ( 6, 15),
          ( 9, 3), ( 9, 6),          ( 9, 12), ( 9, 15),
          (12, 3), (12, 6), (12, 9), (12, 12), (12, 15),
          (14, 4), (14, 14),
          (15, 6), (15, 9), (15, 12)]

FC_pos = [(9, 9)]

# UO₂ fuel assemblies
for n in 1:2
    # corner
    xc = (n - 1) * N * p
    yc = (n - 1) * N * p

    for i in 1:N, j in 1:N
        xo = xc + p / 2 + (i - 1) * p
        yo = yc + p / 2 + (j - 1) * p

        tag = add_circle!(gmsh, (xo, yo), r, lc)
        gmsh.model.geo.addPlaneSurface([tag])

        pos = (i, j)
        if pos in GT_pos
            push!(guidetube, tag)
        elseif pos in FC_pos
            push!(fissionchamber, tag)
        else
            push!(UO₂, tag)
        end
    end
end

MOX_43_pos = Tuple{Int64,Int64}[]
for i in 1:N, j in 1:N
    if i == 1 || i == N || j == 1 || j == N
        push!(MOX_43_pos, (i, j))
    end
end

MOX_87_pos = Tuple{Int64,Int64}[]
for i in 4:N-3, j in 4:N-3
    pos = (i, j)
    if pos in GT_pos
        continue
    elseif pos in FC_pos
        continue
    elseif pos in ((5, 4), (13, 4), (4, 5), (14, 5), (4, 13), (14, 13), (5, 14), (13, 14))
        continue
    else
        push!(MOX_87_pos, pos)
    end
end

MOX_7_pos = Tuple{Int64,Int64}[]
for i in 1:N, j in 1:N
    pos = (i, j)
    if pos in GT_pos
        continue
    elseif pos in FC_pos
        continue
    elseif pos in MOX_43_pos
        continue
    elseif pos in MOX_87_pos
        continue
    # elseif pos in ((5, 4), (13, 4), (4, 5), (14, 5), (4, 13), (14, 13), (5, 14), (13, 14))
    #     push!(MOX_7_pos, pos)
    else
        push!(MOX_7_pos, pos)
    end
end

# MOX fuel assemblies
for n in 1:2
    # corner
    xc = isone(n) ? N * p : 0
    yc = isone(n) ? 0 : N * p

    for i in 1:N, j in 1:N
        xo = xc + p / 2 + (i - 1) * p
        yo = yc + p / 2 + (j - 1) * p

        tag = add_circle!(gmsh, (xo, yo), r, lc)
        gmsh.model.geo.addPlaneSurface([tag])

        pos = (i, j)
        if pos in GT_pos
            push!(guidetube, tag)
        elseif pos in MOX_43_pos
            push!(MOX_43, tag)
        elseif pos in MOX_7_pos
            push!(MOX_7, tag)
        elseif pos in MOX_87_pos
            push!(MOX_87, tag)
        elseif pos in FC_pos
            push!(fissionchamber, tag)
        end
    end
end

s = 3 * N * p
lc = 1.
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

h2oTag = factory.addPlaneSurface(vcat(cl1, guidetube, fissionchamber, UO₂, MOX_43, MOX_7, MOX_87))

pg1 = gmsh.model.geo.addPhysicalGroup(2, guidetube)
pg2 = gmsh.model.geo.addPhysicalGroup(2, fissionchamber)
pg3 = gmsh.model.geo.addPhysicalGroup(2, UO₂)
pg4 = gmsh.model.geo.addPhysicalGroup(2, MOX_43)
pg5 = gmsh.model.geo.addPhysicalGroup(2, MOX_7)
pg6 = gmsh.model.geo.addPhysicalGroup(2, MOX_87)
pg7 = gmsh.model.geo.addPhysicalGroup(2, [h2oTag])
pg8 = gmsh.model.geo.addPhysicalGroup(1, [l1])
pg9 = gmsh.model.geo.addPhysicalGroup(1, [l2])
pg10 = gmsh.model.geo.addPhysicalGroup(1, [l3])
pg11 = gmsh.model.geo.addPhysicalGroup(1, [l4])

gmsh.model.setPhysicalName(2, pg1, "guide-tube")
gmsh.model.setPhysicalName(2, pg2, "fission-chamber")
gmsh.model.setPhysicalName(2, pg3, "UO2")
gmsh.model.setPhysicalName(2, pg4, "MOX_43")
gmsh.model.setPhysicalName(2, pg5, "MOX_7")
gmsh.model.setPhysicalName(2, pg6, "MOX_87")
gmsh.model.setPhysicalName(2, pg7, "water")
gmsh.model.setPhysicalName(1, pg8, "bottom")
gmsh.model.setPhysicalName(1, pg9, "right")
gmsh.model.setPhysicalName(1, pg10, "top")
gmsh.model.setPhysicalName(1, pg11, "left")

# removemos puntos de la geometria duplicados (origenes por ejemplo)
gmsh.model.geo.removeAllDuplicates()

gmsh.model.geo.synchronize()

gmsh.model.mesh.generate(2)

gmsh.write("c5g7.msh")

if !("-nopopup" in ARGS)
    gmsh.fltk.run()
end

gmsh.finalize()

# move to json file format
using Gridap
mshfile = joinpath(@__DIR__,"c5g7.msh")
model = GmshDiscreteModel(mshfile; renumber=true)
Gridap.Io.to_json_file(model, "c5g7.json")
