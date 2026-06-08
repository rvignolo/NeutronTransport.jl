# TODO(feature): support spatially varying cross sections, e.g. via lazy maps over position.
struct CrossSections{
    NGroups,elType,
    T<:Union{Vector{elType},SVector{NGroups,elType}},
    S<:Union{Matrix{elType},SMatrix{NGroups,NGroups,elType}}
}
    name::String
    χ::T
    Σt::T
    νΣf::T
    Σs0::S
    Σs0_sum::T
    fissionable::Bool

    # Candidate data for future formulations:
    # D   :: T1  # diffusion coefficient
    # S   :: T1  # external source
    # Σa  :: T1  # absorption cross section
    # eΣf :: T1  # fission energy release
    # Σs1 :: T2  # first-order scattering matrix
end
const XSs = CrossSections

ngroups(::CrossSections{NGroups}) where {NGroups} = NGroups
eltype(::CrossSections{NGroups,elType}) where {NGroups,elType} = elType
isfissionable(xs::CrossSections) = xs.fissionable

function CrossSections(
    name::AbstractString,
    NGroups::Integer;

    # Future diffusion discretization data.
    # D = nothing,

    # Future fixed-source data.
    # S = nothing,

    Σt = error("Σt has no default, supply it with keyword."),
    Σs0 = error("Σs0 has no default, supply it with keyword."),
    # Σa = nothing, # TODO(feature): derive Σa from Σt and Σs0 when needed.

    # Materials are non-fissionable by default.
    νΣf = nothing,
    χ = nothing,

    # Future anisotropic scattering data.
    # Σs1 = nothing
)
    NGroups > 0 || throw(ArgumentError("`NGroups` must be positive."))

    _check_group_vector(:Σt, Σt, NGroups)
    _check_scattering_matrix(:Σs0, Σs0, NGroups)

    if isnothing(νΣf)
        νΣf = zeros(_float_eltype(eltype(Σt), eltype(Σs0)), NGroups)
    else
        _check_group_vector(:νΣf, νΣf, NGroups)
    end

    if isnothing(χ)
        χ = zeros(_float_eltype(eltype(Σt), eltype(Σs0), eltype(νΣf)), NGroups)
    else
        _check_group_vector(:χ, χ, NGroups)
    end

    promoted = promote_type(eltype(Σt), eltype(Σs0), eltype(νΣf), eltype(χ))
    promoted <: Real || throw(ArgumentError("cross-section data must be real-valued."))
    elType = _float_eltype(promoted)

    use_static = _uses_static_storage(χ, Σt, νΣf, Σs0)

    χ_dense = Vector{elType}(χ)
    Σt_dense = Vector{elType}(Σt)
    νΣf_dense = Vector{elType}(νΣf)
    Σs0_dense = Matrix{elType}(Σs0)
    Σs0_sum_dense = Vector{elType}(undef, NGroups)
    @inbounds for g′ in 1:NGroups
        Σs0_sum_dense[g′] = sum(@view Σs0_dense[g′, :])
    end

    fissionable = any(!iszero, νΣf_dense)
    if fissionable
        sumχ = sum(χ_dense)
        if iszero(sumχ)
            fill!(χ_dense, zero(elType))
            χ_dense[1] = one(elType)
        else
            isapprox(sumχ, one(elType)) ||
                throw(ArgumentError("`χ` must sum to one for fissionable materials."))
            χ_dense ./= sumχ
        end
    else
        fill!(χ_dense, zero(elType))
    end

    χ_out = _cross_section_vector(χ_dense, Val(NGroups), Val(use_static))
    Σt_out = _cross_section_vector(Σt_dense, Val(NGroups), Val(use_static))
    νΣf_out = _cross_section_vector(νΣf_dense, Val(NGroups), Val(use_static))
    Σs0_out = _cross_section_matrix(Σs0_dense, Val(NGroups), Val(use_static))
    Σs0_sum_out = _cross_section_vector(Σs0_sum_dense, Val(NGroups), Val(use_static))

    return CrossSections{NGroups,elType,typeof(χ_out),typeof(Σs0_out)}(
        String(name), χ_out, Σt_out, νΣf_out, Σs0_out, Σs0_sum_out, fissionable
    )
end

function _check_group_vector(name::Symbol, xs, NGroups::Integer)
    xs isa AbstractVector ||
        throw(ArgumentError("`$name` must be a vector with $NGroups entries."))
    length(xs) == NGroups ||
        throw(ArgumentError("`$name` must have length $NGroups; got $(length(xs))."))
    return nothing
end

function _check_scattering_matrix(name::Symbol, xs, NGroups::Integer)
    xs isa AbstractMatrix ||
        throw(ArgumentError("`$name` must be a $NGroups by $NGroups matrix."))
    size(xs) == (NGroups, NGroups) ||
        throw(ArgumentError("`$name` must have size ($NGroups, $NGroups); got $(size(xs))."))
    return nothing
end

_float_eltype(types::Type...) = typeof(float(zero(promote_type(types...))))

_uses_static_storage(xs...) = any(x -> x isa StaticArray, xs)

_cross_section_vector(xs::Vector{T}, ::Val{N}, ::Val{false}) where {N,T} = xs
_cross_section_vector(xs::Vector{T}, ::Val{N}, ::Val{true}) where {N,T} =
    SVector{N,T}(xs)

_cross_section_matrix(xs::Matrix{T}, ::Val{N}, ::Val{false}) where {N,T} = xs
_cross_section_matrix(xs::Matrix{T}, ::Val{N}, ::Val{true}) where {N,T} =
    SMatrix{N,N,T}(xs)
