# IDEA: We can use lazy_map or map to extend each cross section as a function of position
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
    fissionable::Bool

    # D   :: T1
    # S   :: T1
    # Σa  :: T1
    # eΣf :: T1
    # Σs1 :: T2
end
const XSs = CrossSections

ngroups(::CrossSections{NGroups}) where {NGroups} = NGroups
eltype(::CrossSections{NGroups,elType}) where {NGroups,elType} = elType
isfissionable(xs::CrossSections) = xs.fissionable

function CrossSections(
    name::String,
    NGroups::Integer;

    # for future diffusion discretization
    # D = nothing,

    # not yet used
    # S = nothing,

    Σt = error("Σt has no default, supply it with keyword."),
    Σs0 = error("Σs0 has no default, supply it with keyword."),
    # Σa = nothing, # IDEA: can be computed using Σt and Σs0

    # assume not fissionable material
    νΣf = zeros(promote_type(eltype.((Σt, Σs0))...), NGroups),
    χ = zeros(promote_type(eltype.((Σt, νΣf, Σs0))...), NGroups),

    # for future implementations (?)
    # Σs1 = nothing
)

    if iszero(sum(νΣf))
        fissionable = false
        fill!(χ, zero(eltype(χ)))
    else
        fissionable = true
        if iszero(sum(χ))
            χ = zeros(promote_type(eltype.((Σt, νΣf, Σs0))...), NGroups)
            χ[1] = 1
        else
            if !isone(sum(χ))
                error("χ *must* add up to 1.")
            end
        end
    end

    χ, Σt, νΣf = promote(χ, Σt, νΣf)

    if eltype(Σs0) != eltype(χ)
        elType = promote_type(eltype.(χ, Σs0))

        # TODO: convert everything to share element type here y object type (SArray or Array)

        # then, get the types
        T = typeof(χ)
        S = typeof(Σs0)
    else
        elType = eltype(χ)

        # TODO: convert to share object type (SArray or Array) or ar least make sure they match

        T = typeof(χ)
        S = typeof(Σs0)
    end

    return CrossSections{NGroups,elType,T,S}(name, χ, Σt, νΣf, Σs0, fissionable)
end
