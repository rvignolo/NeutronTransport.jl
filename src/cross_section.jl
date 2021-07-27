# IDEA: We can use lazy_map or map to extend each cross section as a function of position
struct CrossSections{G,F,T,T1<:Union{Vector{T},SVector{G,T}},T2<:Union{Matrix{T},SMatrix{G,G,T}}}
    name:: String
    χ   :: T1
    Σt  :: T1
    νΣf :: T1
    Σs0 :: T2

    # D   :: T1
    # S   :: T1
    # Σa  :: T1
    # eΣf :: T1
    # Σs1 :: T2
end

const XSs = CrossSections

ngroups(::CrossSections{G}) where {G} = G
isfissionable(::CrossSections{G,F}) where {G,F} = F
eltype(::CrossSections{G,F,T}) where {G,F,T} = T

function CrossSections(
    name::String,
    G::Integer;

    # for future diffusion discretization
    # D = nothing,

    # not yet used
    # S = nothing,

    Σt = error("Σt has no default, supply it with keyword."),
    Σs0 = error("Σs0 has no default, supply it with keyword."),
    # Σa = nothing, # IDEA: can be computed using Σt and Σs0

    # assume not fissionable material
    νΣf = zeros(promote_type(eltype.((Σt, Σs0))...), G),
    χ = zeros(promote_type(eltype.((Σt, νΣf, Σs0))...), G),

    # for future implementations (?)
    # Σs1 = nothing
)

    if iszero(sum(νΣf))
        F = false
        fill!(χ, zero(eltype(χ)))
    else
        F = true
        if iszero(sum(χ))
            χ = zeros(promote_type(eltype.((Σt, νΣf, Σs0))...), G)
            χ[1] = 1
        else
            if !isone(sum(χ))
                error("χ *must* add up to 1.")
            end
        end
    end

    χ, Σt, νΣf = promote(χ, Σt, νΣf)

    if eltype(Σs0) != eltype(χ)
        T = promote_type(eltype.(χ, Σs0))

        # TODO: convert everything to share element type here y object type (SArray or Array)

        # then, get the types
        T1 = typeof(χ)
        T2 = typeof(Σs0)
    else
        T = eltype(χ)

        # TODO: convert to share object type (SArray or Array) or ar least make sure they match

        T1 = typeof(χ)
        T2 = typeof(Σs0)
    end

    return CrossSections{G,F,T,T1,T2}(name, χ, Σt, νΣf, Σs0)
end
