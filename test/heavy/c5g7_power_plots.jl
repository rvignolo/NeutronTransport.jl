using Printf
using Test
using NeutronTransport

try
    @eval import GridapGmsh: GmshDiscreteModel, gmsh
    @eval using Plots
catch err
    error(
        "C5G7 power plotting requires GridapGmsh.jl and Plots.jl in the active " *
        "Julia load path. Install them in a plotting or benchmark environment " *
        "and rerun this script."
    )
end

const C5G7_RUN_TESTSET = false

include("c5g7_geometry.jl")
include("c5g7.jl")

function parse_output_dir(args)
    for i in eachindex(args)
        if args[i] == "--output-dir" && i < lastindex(args)
            return args[i + 1]
        elseif startswith(args[i], "--output-dir=")
            return split(args[i], "="; limit=2)[2]
        end
    end
    return joinpath(pwd(), "c5g7-power-plots")
end

function c5g7_global_pin_map(values_by_assembly, active_by_assembly)
    n = 2 * C5G7_N_PINS
    values = fill(NaN, n, n)
    active = falses(n, n)

    for assembly in 1:4
        x_offset = isodd(assembly) ? 0 : C5G7_N_PINS
        y_offset = assembly <= 2 ? 0 : C5G7_N_PINS
        for pin_i in 1:C5G7_N_PINS, pin_j in 1:C5G7_N_PINS
            x = x_offset + pin_i
            y = y_offset + pin_j
            values[x, y] = values_by_assembly[assembly, pin_i, pin_j]
            active[x, y] = active_by_assembly[assembly, pin_i, pin_j]
        end
    end

    return values, active
end

function c5g7_power_maps(power_result)
    total_values, active = c5g7_global_pin_map(
        power_result.normalized_pin_powers,
        power_result.active_pin_mask,
    )

    maps = PinPowerMap[
        PinPowerMap(total_values; active=active, title="Total pin power"),
    ]

    for g in 1:size(power_result.normalized_group_pin_powers, 1)
        values, _ = c5g7_global_pin_map(
            @view(power_result.normalized_group_pin_powers[g, :, :, :]),
            power_result.active_pin_mask,
        )
        title = @sprintf("Group %d (%.2f%%)", g, 100 * power_result.group_power_fractions[g])
        push!(maps, PinPowerMap(values; active=active, title=title))
    end

    return maps
end

function add_c5g7_assembly_guides!(plt, subplot)
    split = C5G7_N_PINS + 0.5
    vline!(plt, [split]; subplot=subplot, color=:black, linewidth=1.2, label=false)
    hline!(plt, [split]; subplot=subplot, color=:black, linewidth=1.2, label=false)
    return nothing
end

function c5g7_power_plot(power_result)
    maps = c5g7_power_maps(power_result)
    plt = plot(;
        layout=(2, 4),
        size=(1700, 850),
        plot_title="C5G7 normalized pin-power distribution",
    )

    total_hi = maximum(maps[1].values[maps[1].active])
    plot!(plt, maps[1]; subplot=1, clims=(0, total_hi))
    add_c5g7_assembly_guides!(plt, 1)

    for g in 1:7
        subplot = g + 1
        map = maps[subplot]
        hi = maximum(map.values[map.active])
        plot!(plt, map; subplot=subplot, clims=(0, hi))
        add_c5g7_assembly_guides!(plt, subplot)
    end

    return plt
end

function main(args=ARGS)
    output_dir = parse_output_dir(args)
    mkpath(output_dir)

    result = c5g7_2d_benchmark_result()
    plt = c5g7_power_plot(result.power_result)

    svg_path = joinpath(output_dir, "c5g7_pin_power.svg")
    png_path = joinpath(output_dir, "c5g7_pin_power.png")
    savefig(plt, svg_path)
    savefig(plt, png_path)

    println("wrote ", svg_path)
    println("wrote ", png_path)
    println(@sprintf("keff %.12f", result.sol.keff))
    println(@sprintf(
        "pin power min %.6f max %.6f",
        result.power_result.min_normalized_pin_power,
        result.power_result.max_normalized_pin_power,
    ))

    return (; svg_path, png_path)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
