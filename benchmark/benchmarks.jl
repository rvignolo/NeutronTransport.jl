using BenchmarkTools

include("fixtures.jl")

const SUITE = BenchmarkGroup()

SUITE["constructors"] = BenchmarkGroup()
SUITE["constructors"]["cross_sections"] = @benchmarkable benchmark_single_cell_material()
SUITE["constructors"]["polar_quadrature"] = @benchmarkable TabuchiYamamoto(6)
SUITE["constructors"]["problem"] = @benchmarkable benchmark_single_cell_problem()

SUITE["openmoc"] = BenchmarkGroup()
SUITE["openmoc"]["homogeneous_solve"] = @benchmarkable solve(
    prob; max_iterations=50, max_residual=1e-8
) setup=(prob = benchmark_single_cell_problem(
    bcs=BoundaryConditions(top=Reflective, bottom=Reflective, left=Reflective, right=Reflective),
    material=openmoc_homogeneous_medium(),
))
SUITE["openmoc"]["homogeneous_grid_setup"] = @benchmarkable reflected_problem(
    openmoc_homogeneous_grid_model(),
    [openmoc_homogeneous_medium()]
)
SUITE["openmoc"]["homogeneous_grid_solve_short"] = @benchmarkable solve(
    prob; max_iterations=3, max_residual=0.0
) setup=(prob = reflected_problem(openmoc_homogeneous_grid_model(), [openmoc_homogeneous_medium()]))
SUITE["openmoc"]["reflected_grid_solve_short"] = @benchmarkable solve(
    prob; max_iterations=3, max_residual=0.0
) setup=(prob = reflected_problem(
    openmoc_reflective_grid_model(),
    [openmoc_homogeneous_medium(), openmoc_water_medium()]
))
SUITE["openmoc"]["lattice_grid_solve_short"] = @benchmarkable solve(
    prob; max_iterations=3, max_residual=0.0
) setup=(prob = reflected_problem(
    openmoc_lattice_grid_model(),
    [openmoc_homogeneous_medium(), openmoc_water_medium()]
))
SUITE["openmoc"]["oblong_lattice_solve_short"] = @benchmarkable solve(
    prob; max_iterations=3, max_residual=0.0
) setup=(prob = reflected_problem(
    openmoc_oblong_lattice_grid_model(),
    [openmoc_homogeneous_medium(), openmoc_water_medium(), openmoc_guide_tube_medium()]
))

SUITE["transport"] = BenchmarkGroup()
SUITE["transport"]["optical_length"] = @benchmarkable NeutronTransport.optical_length!(
    prob
) setup=(prob = benchmark_single_cell_problem())
SUITE["transport"]["compute_q"] = @benchmarkable NeutronTransport.compute_q!(
    sol, prob
) setup=begin
    prob = benchmark_single_cell_problem()
    sol = initialized_solution(prob)
end
SUITE["transport"]["compute_phi"] = @benchmarkable NeutronTransport.compute_φ!(
    sol, prob
) setup=begin
    prob = benchmark_single_cell_problem()
    sol = initialized_solution(prob)
end

SUITE["demo"] = BenchmarkGroup()
SUITE["demo"]["pincell_setup"] = @benchmarkable pincell_problem()
SUITE["demo"]["pincell_solve_short"] = @benchmarkable solve(
    prob; max_iterations=3, max_residual=0.0
) setup=(prob = pincell_problem())
SUITE["demo"]["bwr_setup"] = @benchmarkable bwr_problem()
SUITE["demo"]["bwr_solve_short"] = @benchmarkable solve(
    prob; max_iterations=3, max_residual=0.0
) setup=(prob = bwr_problem())
SUITE["demo"]["c5g7_geometry_setup"] = @benchmarkable c5g7_geometry_problem()
SUITE["demo"]["c5g7_geometry_solve_short"] = @benchmarkable solve(
    prob; max_iterations=2, max_residual=0.0
) setup=(prob = c5g7_geometry_problem())

function run_benchmarks(; quick::Bool=false)
    kwargs = quick ? (; samples=1, evals=1, seconds=1) : (; seconds=5)
    return run(SUITE; verbose=true, kwargs...)
end

if abspath(PROGRAM_FILE) == @__FILE__
    results = run_benchmarks(; quick="--quick" in ARGS)
    display(results)
end
