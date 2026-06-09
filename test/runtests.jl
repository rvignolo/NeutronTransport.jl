using NeutronTransport
using StaticArrays
using Test

function run_heavy_tests()
    value = lowercase(get(ENV, "NEUTRONTRANSPORT_RUN_HEAVY_TESTS", "false"))
    return value in ("1", "true", "yes", "on")
end

include("openmoc_references.jl")
include("fixtures.jl")

@testset "NeutronTransport.jl" begin
    include("cross_sections.jl")
    include("polar_quadrature.jl")
    include("openmoc_regressions.jl")
    include("convergence_regressions.jl")
    include("demo_regressions.jl")
    include("plot_recipes.jl")
    include("moc_solver.jl")

    if run_heavy_tests()
        include("heavy/runtests.jl")
    else
        @info(
            "Skipping heavy tests",
            env = "NEUTRONTRANSPORT_RUN_HEAVY_TESTS=true"
        )
    end
end
