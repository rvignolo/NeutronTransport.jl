"""
    fsr_to_cell_values(prob::MoCProblem, fsr_values)

Expand one value per transport flat source region (FSR) to one value per mesh cell.
This is useful for plotting or exporting fields when `cell_to_fsr` coalesces several mesh
cells into the same transport region.
"""
function fsr_to_cell_values(prob::MoCProblem, fsr_values::AbstractVector)
    length(fsr_values) == nregions(prob) ||
        throw(ArgumentError("`fsr_values` must contain one value per transport FSR."))

    cell_values = Vector{eltype(fsr_values)}(undef, length(prob.cell_to_fsr))
    @inbounds for cell in eachindex(prob.cell_to_fsr)
        cell_values[cell] = fsr_values[prob.cell_to_fsr[cell]]
    end

    return cell_values
end

"""
    cell_scalar_flux(sol::MoCSolution, g::Integer)

Return the scalar flux for energy group `g` expanded to one value per mesh cell.
"""
function cell_scalar_flux(sol::MoCSolution, g::Integer)
    return fsr_to_cell_values(sol.prob, sol(g))
end

"""
    CellScalarField(mesh, values; title="")
    CellScalarField(prob::MoCProblem, fsr_values; title="")
    CellScalarField(sol::MoCSolution, g::Integer; title="")

Plot object for a scalar quantity defined on mesh cells. With Plots.jl loaded, use
`plot(CellScalarField(sol, 1))` to render the group-1 scalar flux.
"""
struct CellScalarField{M,V<:AbstractVector}
    mesh::M
    values::V
    title::String
end

function CellScalarField(mesh, values::AbstractVector; title="")
    length(values) == num_cells(mesh) ||
        throw(ArgumentError("`values` must contain one entry per mesh cell."))
    return CellScalarField(mesh, values, string(title))
end

function CellScalarField(prob::MoCProblem, fsr_values::AbstractVector; title="")
    return CellScalarField(
        prob.trackgenerator.mesh, fsr_to_cell_values(prob, fsr_values); title
    )
end

function CellScalarField(sol::MoCSolution, g::Integer; title="")
    title = string(title)
    field_title = isempty(title) ? "Scalar flux group $g" : title
    return CellScalarField(sol.prob.trackgenerator.mesh, cell_scalar_flux(sol, g);
        title=field_title
    )
end

"""
    PinPowerMap(values; active=trues(size(values)), title="")

Plot object for square-lattice pin quantities. `values[i, j]` is interpreted as the value
at pin column `i` and row `j`; inactive pins are hidden from the color scale.
"""
struct PinPowerMap{V<:AbstractMatrix,A<:AbstractMatrix{Bool}}
    values::V
    active::A
    title::String
end

function PinPowerMap(values::AbstractMatrix; active=trues(size(values)), title="")
    size(active) == size(values) ||
        throw(ArgumentError("`active` must have the same size as `values`."))
    return PinPowerMap(values, active, string(title))
end

function _cell_field_coordinates(mesh, values)
    @unpack ordered_cell_nodes, node_coordinates = mesh

    n_cells = length(ordered_cell_nodes)
    max_nodes = maximum(length, ordered_cell_nodes)

    x = Matrix{Float64}(undef, max_nodes + 1, n_cells)
    y = Matrix{Float64}(undef, max_nodes + 1, n_cells)
    fill_z = Vector{Float64}(undef, n_cells)
    fill!(x, NaN)
    fill!(y, NaN)

    @inbounds for (cell, node_ids) in enumerate(ordered_cell_nodes)
        for (j, node_id) in enumerate(node_ids)
            node = node_coordinates[node_id]
            x[j, cell] = node[1]
            y[j, cell] = node[2]
        end

        first_node = node_coordinates[first(node_ids)]
        x[length(node_ids)+1, cell] = first_node[1]
        y[length(node_ids)+1, cell] = first_node[2]
        fill_z[cell] = Float64(values[cell])
    end

    return x, y, fill_z
end

function _masked_pin_values(values, active)
    z = Matrix{Float64}(undef, size(values))
    @inbounds for i in eachindex(values, active)
        z[i] = active[i] ? Float64(values[i]) : NaN
    end
    return z
end

@recipe function plot(field::CellScalarField)
    x, y, z = _cell_field_coordinates(field.mesh, field.values)

    seriestype := :shape
    fill_z --> z
    linecolor --> :transparent
    linewidth --> 0
    legend --> false
    colorbar --> true
    aspect_ratio --> :equal
    framestyle --> :box
    seriescolor --> :viridis
    title --> field.title

    return x, y
end

@recipe function plot(map::PinPowerMap)
    z = _masked_pin_values(map.values, map.active)

    seriestype := :heatmap
    legend --> false
    colorbar --> true
    aspect_ratio --> :equal
    framestyle --> :box
    seriescolor --> :viridis
    ticks --> false
    title --> map.title

    return axes(z, 1), axes(z, 2), permutedims(z)
end
