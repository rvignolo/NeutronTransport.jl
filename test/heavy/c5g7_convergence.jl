using Printf
using Test
using NeutronTransport

import GridapGmsh: GmshDiscreteModel, gmsh

const C5G7_RUN_TESTSET = false

include("c5g7_geometry.jl")
include("c5g7.jl")

const C5G7_CONVERGENCE_CASES = (
    (
        name="baseline_ty",
        pin_lc=0.25,
        polar_quadrature=TabuchiYamamoto,
        max_residual=1e-4,
    ),
    (
        name="baseline_ty_residual_1e-5",
        pin_lc=0.25,
        polar_quadrature=TabuchiYamamoto,
        max_residual=1e-5,
        max_iterations=1000,
    ),
    (
        name="equal_angle",
        pin_lc=0.25,
        polar_quadrature=EqualAngle,
        max_residual=1e-4,
    ),
    (
        name="equal_angle_residual_1e-5",
        pin_lc=0.25,
        polar_quadrature=EqualAngle,
        max_residual=1e-5,
        max_iterations=1000,
    ),
    (
        name="baseline_ty_pin_lc_0_15",
        pin_lc=0.15,
        polar_quadrature=TabuchiYamamoto,
        max_residual=1e-4,
    ),
    (
        name="baseline_ty_pin_lc_0_10",
        pin_lc=0.10,
        polar_quadrature=TabuchiYamamoto,
        max_residual=1e-4,
        segmentize_kwargs=(k=20, rtol=2 * sqrt(eps(Float64))),
    ),
    (
        name="equal_angle_pin_lc_0_15",
        pin_lc=0.15,
        polar_quadrature=EqualAngle,
        max_residual=1e-4,
    ),
)

function c5g7_case_kwargs(case)
    names = Tuple(k for k in keys(case) if k != :name)
    values = Tuple(getfield(case, k) for k in names)
    return NamedTuple{names}(values)
end

function c5g7_case_by_name(name)
    for case in C5G7_CONVERGENCE_CASES
        case.name == name && return case
    end
    throw(ArgumentError("unknown C5G7 convergence case: $name"))
end

function parse_case_names(args)
    names = String[]
    for arg in args
        if startswith(arg, "--case=")
            push!(names, split(arg, "="; limit=2)[2])
        elseif startswith(arg, "--cases=")
            append!(names, split(split(arg, "="; limit=2)[2], ","))
        end
    end

    return isempty(names) ? [case.name for case in C5G7_CONVERGENCE_CASES] : names
end

function c5g7_convergence_result(case)
    reference = C5G7_2D_REFERENCE
    result = c5g7_2d_benchmark_result(; c5g7_case_kwargs(case)...)
    power_result = result.power_result

    material_volume_relerr = abs.(
        (result.material_volumes .- result.reference_material_volumes) ./
        result.reference_material_volumes
    )
    assembly_power_relerr = abs.(
        (collect(values(power_result.assembly_powers)) .-
         collect(values(reference.assembly_powers))) ./
        collect(values(reference.assembly_powers))
    )
    max_pin_power_relerr = abs(
        power_result.max_normalized_pin_power - reference.max_normalized_pin_power
    ) / reference.max_normalized_pin_power
    min_pin_power_relerr = abs(
        power_result.min_normalized_pin_power - reference.min_normalized_pin_power
    ) / reference.min_normalized_pin_power

    return (;
        name=case.name,
        pin_lc=result.settings.pin_lc,
        polar_quadrature=string(result.settings.polar_quadrature),
        max_residual=result.settings.max_residual,
        keff=result.sol.keff,
        keff_error_pcm=result.keff_error_pcm,
        iterations=result.sol.iterations,
        residual=result.sol.residual,
        max_material_volume_relerr=maximum(material_volume_relerr),
        max_pin_power_relerr,
        min_pin_power_relerr,
        max_assembly_power_relerr=maximum(assembly_power_relerr),
        n_cells=result.n_cells,
        n_segments=result.n_segments,
        timings=result.timings,
    )
end

function print_c5g7_convergence_result(result)
    @printf(
        "%-28s pin_lc=%4.2f polar=%-15s residual=%7.1e keff=%.9f pcm=%+8.1f iters=%4d cells=%7d segments=%8d solve=%7.1fs vol=%6.3f%% pin_max=%6.3f%% pin_min=%6.3f%% asm=%6.3f%%\n",
        result.name,
        result.pin_lc,
        result.polar_quadrature,
        result.max_residual,
        result.keff,
        result.keff_error_pcm,
        result.iterations,
        result.n_cells,
        result.n_segments,
        result.timings.solve,
        100 * result.max_material_volume_relerr,
        100 * result.max_pin_power_relerr,
        100 * result.min_pin_power_relerr,
        100 * result.max_assembly_power_relerr,
    )
end

function main(args=ARGS)
    results = []
    for name in parse_case_names(args)
        case = c5g7_case_by_name(name)
        result = c5g7_convergence_result(case)
        print_c5g7_convergence_result(result)
        push!(results, result)
    end
    return results
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
