using NeutronTransport
import Gridap: DiscreteModelFromFile

jsonfile = joinpath(@__DIR__,"subcrit-assy.json")
geometry = DiscreteModelFromFile(jsonfile)

# number of azimuthal angles
nφ = 4

# azimuthal spacing
δ = 0.1

# boundary conditions
bcs = BoundaryConditions(top=Reflective, bottom=Reflective, left=Reflective, right=Reflective)

# initialize track generator
tg = TrackGenerator(geometry, nφ, δ, bcs=bcs)

# perform ray tracing
trace!(tg)

# proceed to segmentation
segmentize!(tg)

# polar quadrature
pq = TabuchiYamamoto(6)

# materials
water = CrossSections("water", 7;
    νΣf = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    Σt  = [1.592060E-01, 4.129700E-01, 5.903100E-01, 5.843500E-01, 7.180000E-01, 1.254450E+00, 2.650380E+00],
    Σs0 = [4.447770E-02 1.134000E-01 7.234700E-04 3.749900E-06 5.318400E-08 0. 0.; 
            0. 2.823340E-01 1.299400E-01 6.234000E-04 4.800200E-05 7.448600E-06 1.045500E-06; 
            0. 0. 3.452560E-01 2.245700E-01 1.699900E-02 2.644300E-03 5.034400E-04;
            0. 0. 0. 9.102840E-02 4.155100E-01 6.373200E-02 1.213900E-02;
            0. 0. 0. 7.143700E-05 1.391380E-01 5.118200E-01 6.122900E-02; 
            0. 0. 0. 0. 2.215700E-03 6.999130E-01 5.373200E-01;
            0. 0. 0. 0. 0. 1.324400E-01 2.480700E+00]
)

source = CrossSections("source", 7;
    νΣf = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    Σt  = [1.592060E-01, 4.129700E-01, 5.903100E-01, 5.843500E-01, 7.180000E-01, 1.254450E+00, 2.650380E+00],
    Σs0 = [4.447770E-02 1.134000E-01 7.234700E-04 3.749900E-06 5.318400E-08 0. 0.; 
            0. 2.823340E-01 1.299400E-01 6.234000E-04 4.800200E-05 7.448600E-06 1.045500E-06; 
            0. 0. 3.452560E-01 2.245700E-01 1.699900E-02 2.644300E-03 5.034400E-04;
            0. 0. 0. 9.102840E-02 4.155100E-01 6.373200E-02 1.213900E-02;
            0. 0. 0. 7.143700E-05 1.391380E-01 5.118200E-01 6.122900E-02; 
            0. 0. 0. 0. 2.215700E-03 6.999130E-01 5.373200E-01;
            0. 0. 0. 0. 0. 1.324400E-01 2.480700E+00]
)


fuel = CrossSections("fuel", 7;
νΣf = [2.005998E-02, 2.027303E-03, 1.570599E-02, 4.518301E-02, 4.334208E-02, 2.020901E-01, 5.257105E-01],
Σt  = [1.779490E-01, 3.298050E-01, 4.803880E-01, 5.543670E-01, 3.118010E-01, 3.951680E-01, 5.644060E-01],
Σs0 = [1.275370E-01 4.237800E-02 9.437400E-06 5.516300E-09 0. 0. 0.;
        0. 3.244560E-01 1.631400E-03 3.142700E-09 0. 0. 0.;
        0. 0. 4.509400E-01 2.679200E-03 0. 0. 0.;
        0. 0. 0. 4.525650E-01 5.566400E-03 0. 0.; 
        0. 0. 0. 1.252500E-04 2.714010E-01 1.025500E-02 1.002100E-08;
        0. 0. 0. 0. 1.296800E-03 2.658020E-01 1.680900E-02;
        0. 0. 0. 0. 0. 8.545800E-03 2.730800E-01]
)

xs = [water, source, fuel]

# define the problem
prob = MoCProblem(tg, pq, xs)

# define fixed source material
fixed_sources = set_fixed_source_material(prob, "source", 3, 1e-2)

# solve
sol = solve(prob, fixed_sources, debug=true, max_iterations=2000)

import Gridap: writevtk
import Gridap.Geometry: get_triangulation
trian = get_triangulation(tg.mesh.model)
writevtk(trian, "subcrit-assy", cellfields=[string(g) => sol(g) for g in 1:7])