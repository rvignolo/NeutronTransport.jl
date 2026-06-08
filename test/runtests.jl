using NeutronTransport
using StaticArrays
using Test

include("openmoc_references.jl")
include("fixtures.jl")

@testset "NeutronTransport.jl" begin
    include("cross_sections.jl")
    include("polar_quadrature.jl")
    include("openmoc_regressions.jl")
    include("demo_regressions.jl")
    include("moc_solver.jl")
end
