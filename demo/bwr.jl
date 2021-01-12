using NeutronTransport
using Gridap

jsonfile = joinpath(@__DIR__,"bwr.json")
geometry = DiscreteModelFromFile(jsonfile)

# number of azimuthal angles
nφ = 16

# azimuthal spacing
δ = 0.01

# boundary conditions
bcs = BoundaryConditions(top=Reflective, bottom=Reflective, left=Reflective, right=Reflective)
# bcs = BoundaryConditions(top=Vaccum, bottom=Vaccum, left=Vaccum, right=Vaccum)
# bcs = BoundaryConditions(top=Periodic, bottom=Periodic, left=Periodic, right=Periodic)

# initialize track generator
tg = TrackGenerator(geometry, nφ, δ, bcs=bcs)

# perform ray tracing
trace!(tg)

# proceed to segmentation
segmentize!(tg)

# polar quadrature
pq = TabuchiYamamoto(6)

# materials
pin = CrossSections(2;
    νΣf = [1.86278e-2, 3.44137e-1],
    Σt  = [3.62022e-1, 5.72155e-1],
    Σs0 = [3.33748e-1  6.64881e-4; 0.0e-0 3.80898e-1]
)

cladding = CrossSections(2;
    Σt  = [2.74144e-1, 2.80890e-1],
    Σs0 = [2.72377e-1 1.90838e-4; 0.0e-0 2.77230e-1]
)

water = CrossSections(2;
    Σt  = [6.40711e-1, 1.69131e-0],
    Σs0 = [6.07382e-1 3.31316e-2; 0.0e-0 1.68428e-0]
)

pin_gd = CrossSections(2;
    νΣf = [1.79336e-2, 1.57929e-1],
    Σt  = [3.71785e-1, 1.75000e-0],
    Σs0 = [3.38096e-1  6.92807e-4; 0.0e-0 3.83204e-1]
)

materials = Dict("pin" => pin, "cladding" => cladding, "water" => water, "pin-gd" => pin_gd)

# define the problem
prob = MoCProblem(tg, pq, materials)

# solve
sol = solve(prob)

φ1 = Vector{Float64}(undef, Int64(length(sol.φ)/2));
j = 0
for i in 1:length(sol.φ)
    if isodd(i)
        j += 1
        φ1[j] = sol.φ[i]
    end
end

φ2 = Vector{Float64}(undef, Int64(length(sol.φ)/2));
j = 0
for i in 1:length(sol.φ)
    if iseven(i)
        j += 1
        φ2[j] = sol.φ[i]
    end
end

trian = Gridap.Geometry.get_triangulation(tg.mesh.model)

writevtk(trian, "fluxes", cellfields=["grupo2" => φ2, "grupo1" => φ1])