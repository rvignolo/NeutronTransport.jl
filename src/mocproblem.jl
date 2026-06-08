import RayTracing: num_dims, num_cells

"""
    MoCProblem{Dim,NRegions,NGroups,elType,G<:TrackGenerator,Q<:Quadrature,Xs<:XSs} <: TransportProblem{Dim,NRegions,NGroups}

Method of Characteristics problem built from a ray-traced mesh, polar quadrature, and
material cross sections.
"""
struct MoCProblem{Dim,NRegions,NGroups,elType,G<:TrackGenerator,Q<:Quadrature,Xs<:XSs} <: TransportProblem{Dim,NRegions,NGroups}
    nφ::Int
    nψ::Int

    trackgenerator::G
    quadrature::Q

    xss::Vector{Xs}
    cell_to_fsr::Vector{Int32}  # Mesh cell id -> transport FSR id.
    fsr_tag::Vector{Int32}      # Transport FSR id -> cross-section index in `xss`.
    volumes::Vector{elType}     # Transport FSR volumes.
    optical_lengths::Vector{Matrix{elType}}
    attenuation_factors::Vector{Array{elType,3}}

    function MoCProblem{Dim,NRegions,NGroups,elType}(
        nφ, nψ, tg::G, quad::Q, xss::Vector{X}, cell_to_fsr, fsr_tag, volumes,
        optical_lengths, attenuation_factors
    ) where {Dim,NRegions,NGroups,elType,G,Q,X}
        return new{Dim,NRegions,NGroups,elType,G,Q,X}(
            nφ, nψ, tg, quad, xss, cell_to_fsr, fsr_tag, volumes,
            optical_lengths, attenuation_factors
        )
    end
end

function MoCProblem(
    tg::TrackGenerator,
    polar_quadrature::PolarQuadrature,
    xss::Vector{<:XSs};
    cell_to_fsr=nothing,
    fsr_to_xs=nothing
)
    @unpack mesh, azimuthal_quadrature, n_total_tracks = tg

    Dim = num_dims(mesh)

    NGroups = ngroups(first(xss))
    if !all(g -> isequal(g, NGroups), ngroups.(xss))
        error("all `CrossSections` *must* have the same number of energy groups.")
    end
    elType = promote_type(eltype.(xss)...)

    n_cells = num_cells(mesh)
    cell_xs = cell_material_ids(tg, xss)
    cell_to_fsr, NRegions = _resolve_cell_to_fsr(cell_to_fsr, n_cells)
    fsr_tag = _resolve_fsr_to_xs(fsr_to_xs, cell_to_fsr, cell_xs, NRegions, length(xss))
    volumes = _aggregate_fsr_volumes(tg.volumes, cell_to_fsr, NRegions, elType)

    n_polar_half_count = n_polar_half(polar_quadrature)

    nφ = NRegions * NGroups
    nψ = n_total_tracks * 2 * n_polar_half_count * NGroups  # forward/backward directions.

    quadrature = Quadrature(azimuthal_quadrature, polar_quadrature)

    optical_lengths = [Matrix{elType}(undef, NGroups, 0) for _ in 1:n_total_tracks]
    attenuation_factors = [
        Array{elType,3}(undef, NGroups, n_polar_half_count, 0) for _ in 1:n_total_tracks
    ]

    return MoCProblem{Dim,NRegions,NGroups,elType}(
        nφ, nψ, tg, quadrature, xss, cell_to_fsr, fsr_tag, volumes,
        optical_lengths, attenuation_factors
    )
end

"""
    cell_material_ids(tg::TrackGenerator, xss::Vector{<:CrossSections})

Return one cross-section index per mesh cell by matching Gridap material labels against the
names in `xss`.
"""
function cell_material_ids(tg::TrackGenerator, xss::Vector{<:XSs})
    Dim = num_dims(tg.mesh)
    face_labeling = get_face_labeling(tg.mesh.model)
    tag_to_idx = Dict{Int,Int32}()
    for (i, xs) in enumerate(xss)
        # Resolve each material name to the corresponding Gridap cell tag.
        gridap_tag = get_tag_from_name(face_labeling, xs.name)

        # Map Gridap tags to the NeutronTransport cross-section index.
        push!(tag_to_idx, Int(gridap_tag) => Int32(i))
    end

    xs_names = [xs.name for xs in xss]
    cell_tag = get_face_tag(face_labeling, xs_names, Dim)
    cell_xs = Vector{Int32}(undef, length(cell_tag))
    for (i, gridap_tag) in enumerate(cell_tag)
        cell_xs[i] = tag_to_idx[Int(gridap_tag)]
    end

    return cell_xs
end

function _resolve_cell_to_fsr(cell_to_fsr, n_cells::Integer)
    if isnothing(cell_to_fsr)
        return Int32.(1:n_cells), n_cells
    end

    cell_to_fsr isa AbstractVector{<:Integer} ||
        throw(ArgumentError("`cell_to_fsr` must be an integer vector."))
    length(cell_to_fsr) == n_cells ||
        throw(ArgumentError("`cell_to_fsr` must have one entry per mesh cell."))

    mapping = Int32.(cell_to_fsr)
    all(>=(1), mapping) ||
        throw(ArgumentError("`cell_to_fsr` ids must be positive 1-based integers."))

    n_fsrs = Int(maximum(mapping))
    seen = falses(n_fsrs)
    for fsr in mapping
        seen[fsr] = true
    end
    all(seen) ||
        throw(ArgumentError("`cell_to_fsr` ids must be contiguous from 1 to N."))

    return mapping, n_fsrs
end

function _resolve_fsr_to_xs(fsr_to_xs, cell_to_fsr, cell_xs, n_fsrs::Integer, n_xs::Integer)
    if !isnothing(fsr_to_xs)
        fsr_to_xs isa AbstractVector{<:Integer} ||
            throw(ArgumentError("`fsr_to_xs` must be an integer vector."))
        length(fsr_to_xs) == n_fsrs ||
            throw(ArgumentError("`fsr_to_xs` must have one entry per FSR."))

        mapping = Int32.(fsr_to_xs)
        all(xs_idx -> 1 <= xs_idx <= n_xs, mapping) ||
            throw(ArgumentError("`fsr_to_xs` entries must index `xss`."))
        return mapping
    end

    mapping = zeros(Int32, n_fsrs)
    for cell in eachindex(cell_to_fsr)
        fsr = cell_to_fsr[cell]
        xs_idx = cell_xs[cell]

        if iszero(mapping[fsr])
            mapping[fsr] = xs_idx
        elseif mapping[fsr] != xs_idx
            throw(ArgumentError(
                "FSR $fsr contains cells with multiple materials. Pass explicit " *
                "`fsr_to_xs` only if you intentionally use homogenized cross sections."
            ))
        end
    end

    return mapping
end

function _aggregate_fsr_volumes(cell_volumes, cell_to_fsr, n_fsrs::Integer, ::Type{T}) where {T}
    length(cell_volumes) == length(cell_to_fsr) ||
        throw(ArgumentError("track-generated volumes must have one entry per mesh cell."))

    volumes = zeros(T, n_fsrs)
    for cell in eachindex(cell_to_fsr)
        volumes[cell_to_fsr[cell]] += cell_volumes[cell]
    end

    return volumes
end

dimension(::MoCProblem{Dim}) where {Dim} = Dim
nregions(::MoCProblem{Dim,NRegions}) where {Dim,NRegions} = NRegions
ngroups(::MoCProblem{Dim,NRegions,NGroups}) where {Dim,NRegions,NGroups} = NGroups
eltype(::MoCProblem{Dim,NRegions,NGroups,elType}) where {Dim,NRegions,NGroups,elType} = elType

function show(io::IO, prob::MoCProblem)
    println(io, "  Problem Dimension: ", dimension(prob))
    println(io, "  Number of regions: ", nregions(prob))
    print(io, "  Number of energy groups: ", ngroups(prob))
end

@inline fsr_id(prob::MoCProblem, cell_idx::Integer) = prob.cell_to_fsr[cell_idx]

@inline function getxs(prob::MoCProblem, region_idx::Integer)
    @unpack xss, fsr_tag = prob
    @inbounds xs_idx = getindex(fsr_tag, region_idx)
    @inbounds xs = getindex(xss, xs_idx)
    return xs
end

# Scalar flux and source vectors are region-major: all energy groups for one region are
# contiguous. A group-major layout would use `(g - 1) * NRegions + i`.
macro region_index(i, g)
    ex = quote
        ($i - 1) * NGroups + $g
    end
    return esc(ex)
end

# These index helpers are macros because they expand against local loop constants such as
# `NGroups` and `n_polar_half_count`.

macro angular_index(t, d, p, g)
    ex = quote
        ($t - 1) * 2 * n_polar_half_count * NGroups +
        $d * n_polar_half_count * NGroups + ($p - 1) * NGroups + $g
    end
    return esc(ex)
end

macro reduced_angular_index(p, g)
    ex = quote
        ($p - 1) * NGroups + $g
    end
    return esc(ex)
end
