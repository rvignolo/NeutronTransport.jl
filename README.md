# NeutronTransport

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://rvignolo.github.io/NeutronTransport.jl/dev/)
[![Build Status](https://github.com/rvignolo/NeutronTransport.jl/workflows/CI/badge.svg)](https://github.com/rvignolo/NeutronTransport.jl/actions)

**NeutronTransport** is a reactor physics code that solves the steady-state multigroup neutron transport equation by means of the Method of Characteristics approximation over unstructured grids.

## Examples

Let's solve the transport equation over a simple *pin-cell* geometry:

```julia
using NeutronTransport
using Gridap

# load mesh from json file
jsonfile = joinpath(@__DIR__, "demo", "pincell.json")
geometry = DiscreteModelFromFile(jsonfile)

# number of azimuthal angles
nφ = 32

# azimuthal spacing
δ = 1e-2

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

# materials container
materials = Dict("pin" => pin, "cladding" => cladding, "water" => water)

# define the problem
prob = MoCProblem(tg, pq, materials)

# solve
sol = solve(prob)
```

| ![](https://github.com/rvignolo/NeutronTransport.jl/blob/main/demo/pincell-fast.png)  |  ![](https://github.com/rvignolo/NeutronTransport.jl/blob/main/demo/pincell-thermal.png) |
|:-------------:|:-------------:|
| Fast flux | Thermal flux |

Now, let's solve the transport equation in a 4 by 4 BWR lattice with two Gadolinium pins:
