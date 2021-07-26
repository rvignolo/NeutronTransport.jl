# IDEA: We can use lazy_map or map to extend each cross section as a function of position
struct CrossSections{G,T,T1<:Union{Vector{T},SVector{G,T}},T2<:Union{Matrix{T},SMatrix{G,G,T}}}
    name:: String
    D   :: T1
    S   :: T1
    χ   :: T1
    Σt  :: T1
    Σa  :: T1
    νΣf :: T1
    eΣf :: T1
    Σs0 :: T2
    Σs1 :: T2
end

const XSs = CrossSections

ngroups(::CrossSections{NGroups}) where {NGroups} = NGroups
eltype(::CrossSections{NGroups,T}) where {NGroups,T} = T

function CrossSections(
    name::String,
    G::Integer;
    D=zeros(G), S=zeros(G), χ=zeros(G),
    Σt=zeros(G), Σa=zeros(G), νΣf=zeros(G), eΣf=zeros(G), Σs0=zeros(G, G), Σs1=zeros(G, G)
)
    # TODO: check que se cumplan las relaciones necesarias, ej: Σt no es indep. de Σa y Σs
    # TODO: check que se provean las necesarias si o si.
    # TODO: si hay al menos un SArray, convertir todos a ese tipo.
    # TODO: check correct sizes of all objects
    T = eltype(Σt)
    T1 = typeof(Σt)
    T2 = typeof(Σs0)
    if iszero(sum(χ))
        χ[begin] = one(eltype(χ))
    else
        if !isone(sum(χ))
            error("χ does not add up to 1.")
        end
    end
    return CrossSections{G,T,T1,T2}(name, D, S, χ, Σt, Σa, νΣf, eΣf, Σs0, Σs1)
end