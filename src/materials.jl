# IDEA: We can use lazy_map or map to extend each cross section as a function of position

# IDEA: material name could be inside the `CrossSections` struct, right?
struct CrossSections{
    G,T1<:Union{Vector{T},SVector{G,T}} where T,T2<:Union{Matrix{T},SMatrix{G,G,T}} where T
}
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

ngroups(::CrossSections{G}) where {G} = G

function CrossSections(
    g::Integer;
    D=zeros(g), S=zeros(g), χ=zeros(g),
    Σt=zeros(g), Σa=zeros(g), νΣf=zeros(g), eΣf=zeros(g), Σs0=zeros(g, g), Σs1=zeros(g, g)
)
    # TODO: check que se cumplan las relaciones necesarias, ej: Σt no es indep. de Σa y Σs
    # TODO: check que se provean las necesarias si o si.
    # TODO: si hay al menos un SArray, convertir todos a ese tipo.
    # TODO: check correct sizes of all objects
    T1 = typeof(Σt)
    T2 = typeof(Σs0)
    if iszero(sum(χ))
        χ[begin] = one(eltype(χ))
    else
        if !isone(sum(χ))
            error("χ does not add up to 1.")
        end
    end
    return CrossSections{g,T1,T2}(D, S, χ, Σt, Σa, νΣf, eΣf, Σs0, Σs1)
end