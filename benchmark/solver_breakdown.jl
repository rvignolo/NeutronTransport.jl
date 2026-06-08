using Printf
using Statistics

include("fixtures.jl")

const DEFAULT_CASE = "c5g7-demo"
const DEFAULT_SAMPLES = 10
const BENCHMARK_SINK = Ref{Any}()

function arg_value(args, name, default)
    prefix = name * "="
    for i in eachindex(args)
        arg = args[i]
        if startswith(arg, prefix)
            return split(arg, "="; limit=2)[2]
        elseif arg == name && i < lastindex(args)
            return args[i + 1]
        end
    end
    return default
end

has_flag(args, name) = name in args

function make_problem(case::AbstractString)
    if case == "c5g7-demo"
        return c5g7_geometry_problem(; n_azim=8, spacing=0.25, n_polar=4)
    elseif case == "bwr"
        return bwr_problem(; n_azim=8, spacing=0.12, n_polar=4)
    elseif case == "openmoc-grid"
        return reflected_problem(
            openmoc_homogeneous_grid_model(),
            [openmoc_homogeneous_medium()];
            n_azim=16,
            spacing=0.125,
            n_polar=4,
        )
    else
        error("unknown benchmark case `$case`; use c5g7-demo, bwr, or openmoc-grid")
    end
end

function prepared_solution(prob; parallel::Bool=false)
    sol = NeutronTransport.MoCSolution{eltype(prob)}(prob)
    NeutronTransport.optical_length!(prob)
    NeutronTransport.set_uniform_φ!(sol, one(eltype(prob)))
    NeutronTransport.set_uniform_Q!(sol, zero(eltype(prob)))
    NeutronTransport.set_uniform_start_boundary_ψ!(sol, zero(eltype(prob)))
    NeutronTransport.update_boundary_ψ!(sol)
    NeutronTransport.normalize_fluxes!(sol, prob; parallel=parallel)
    NeutronTransport.compute_q!(sol, prob; parallel=parallel)
    return sol
end

function prepared_residual_solution(prob; parallel::Bool=false)
    sol = prepared_solution(prob; parallel=parallel)
    NeutronTransport.compute_φ!(sol, prob; parallel=parallel)
    NeutronTransport.residual(sol, prob; parallel=parallel)
    return sol
end

function count_segments(prob)
    n_segments = 0
    for track in prob.trackgenerator.tracks_by_uid
        n_segments += length(track.segments)
    end
    return n_segments
end

function format_time(ns::Real)
    if ns < 1_000
        return @sprintf("%.2f ns", ns)
    elseif ns < 1_000_000
        return @sprintf("%.2f us", ns / 1_000)
    elseif ns < 1_000_000_000
        return @sprintf("%.2f ms", ns / 1_000_000)
    else
        return @sprintf("%.2f s", ns / 1_000_000_000)
    end
end

function measure!(kernel, results, name, samples::Integer, setup)
    times = Vector{Float64}(undef, samples)
    for sample in 1:samples
        GC.gc()
        state = setup()
        t0 = time_ns()
        BENCHMARK_SINK[] = kernel(state)
        times[sample] = time_ns() - t0
    end
    t = median(times)
    results[name] = t
    println(rpad(name, 32), format_time(t))
    return nothing
end

function run_case(case; samples::Int=DEFAULT_SAMPLES)
    results = Dict{String,Float64}()

    prob = make_problem(case)
    println("case: ", case)
    println("samples: ", samples)
    println("threads: ", Base.Threads.nthreads())
    println("regions: ", NeutronTransport.nregions(prob))
    println("groups: ", NeutronTransport.ngroups(prob))
    println("tracks: ", length(prob.trackgenerator.tracks_by_uid))
    println("segments: ", count_segments(prob))
    println()

    measure!(results, "problem_setup", samples, () -> nothing) do _
        make_problem(case)
    end

    measure!(results, "optical_length!", samples, () -> make_problem(case)) do prob
        NeutronTransport.optical_length!(prob)
    end

    measure!(results, "normalize_fluxes!", samples, () -> begin
        prob = make_problem(case)
        sol = prepared_solution(prob)
        return (; prob, sol)
    end) do state
        NeutronTransport.normalize_fluxes!(state.sol, state.prob)
    end

    measure!(results, "total_fission_source", samples, () -> begin
        prob = make_problem(case)
        sol = prepared_solution(prob)
        return (; prob, sol)
    end) do state
        NeutronTransport.total_fission_source(state.sol, state.prob)
    end

    measure!(results, "compute_q!", samples, () -> begin
        prob = make_problem(case)
        sol = prepared_solution(prob)
        return (; prob, sol)
    end) do state
        NeutronTransport.compute_q!(state.sol, state.prob)
    end

    measure!(results, "compute_φ!", samples, () -> begin
        prob = make_problem(case)
        sol = prepared_solution(prob)
        return (; prob, sol)
    end) do state
        NeutronTransport.compute_φ!(state.sol, state.prob)
    end

    measure!(results, "residual", samples, () -> begin
        prob = make_problem(case)
        sol = prepared_residual_solution(prob)
        return (; prob, sol)
    end) do state
        NeutronTransport.residual(state.sol, state.prob)
    end

    measure!(results, "solve_3_iter", samples, () -> make_problem(case)) do prob
        solve(prob; max_iterations=3, max_residual=0.0)
    end

    if Base.Threads.nthreads() > 1
        measure!(results, "total_fission_source_parallel", samples, () -> begin
            prob = make_problem(case)
            sol = prepared_solution(prob; parallel=true)
            return (; prob, sol)
        end) do state
            NeutronTransport.total_fission_source(state.sol, state.prob; parallel=true)
        end

        measure!(results, "compute_q!_parallel", samples, () -> begin
            prob = make_problem(case)
            sol = prepared_solution(prob; parallel=true)
            return (; prob, sol)
        end) do state
            NeutronTransport.compute_q!(state.sol, state.prob; parallel=true)
        end

        measure!(results, "compute_φ!_parallel", samples, () -> begin
            prob = make_problem(case)
            sol = prepared_solution(prob; parallel=true)
            return (; prob, sol)
        end) do state
            NeutronTransport.compute_φ!(state.sol, state.prob; parallel=true)
        end

        measure!(results, "residual_parallel", samples, () -> begin
            prob = make_problem(case)
            sol = prepared_residual_solution(prob; parallel=true)
            return (; prob, sol)
        end) do state
            NeutronTransport.residual(state.sol, state.prob; parallel=true)
        end

        measure!(results, "solve_3_iter_parallel", samples, () -> make_problem(case)) do prob
            solve(prob; max_iterations=3, max_residual=0.0, parallel=true)
        end
    end

    println()
    print_breakdown(results)
    print_parallel_ratios(results)

    return results
end

function print_breakdown(results)
    names = [
        "normalize_fluxes!",
        "compute_q!",
        "compute_φ!",
        "total_fission_source",
        "residual",
    ]
    total = sum(results[name] for name in names)
    println("serial one-iteration kernel share:")
    for name in names
        pct = 100 * results[name] / total
        println("  ", rpad(name, 24), @sprintf("%6.2f%%", pct))
    end
end

function print_parallel_ratios(results)
    pairs = [
        ("total_fission_source", "total_fission_source_parallel"),
        ("compute_q!", "compute_q!_parallel"),
        ("compute_φ!", "compute_φ!_parallel"),
        ("residual", "residual_parallel"),
        ("solve_3_iter", "solve_3_iter_parallel"),
    ]
    ratios = filter(p -> haskey(results, p[2]), pairs)
    isempty(ratios) && return nothing

    println()
    println("parallel speedups:")
    for (serial_name, parallel_name) in ratios
        speedup = results[serial_name] / results[parallel_name]
        println("  ", rpad(serial_name, 24), @sprintf("%6.2fx", speedup))
    end
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    args = copy(ARGS)
    samples = parse(Int, arg_value(args, "--samples", string(DEFAULT_SAMPLES)))
    case = arg_value(args, "--case", DEFAULT_CASE)
    if has_flag(args, "--quick")
        samples = min(samples, 3)
    end
    run_case(case; samples=samples)
end
