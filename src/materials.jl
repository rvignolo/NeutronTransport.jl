# IDEA: We can use lazy_map or map to extend each cross section as a function of position

struct CrossSections{G,T1<:Union{Vector{T},SVector{G,T}} where T,
    T2<:Union{Matrix{T},SMatrix{G,G,T}} where T}
    D::T1
    S::T1
    χ::T1
    Σt::T1
    Σa::T1
    νΣf::T1
    eΣf::T1
    Σs0::T2
    Σs1::T2
end

ngroups(::CrossSections{G}) where {G} = G

function CrossSections(g::Int64,
    T=Float64;
    D=zeros(T, g), S=zeros(T, g), χ=nothing, Σt=zeros(T, g), Σa=zeros(T, g),
    νΣf=zeros(T, g), eΣf=zeros(T, g), Σs0=zeros(T, g, g), Σs1=zeros(T, g, g)
)
    # TODO: check que se cumplan las relaciones necesarias, ej: Σt no es indep. de Σa y Σs
    # TODO: check que se provean las necesarias si o si.
    # TODO: si hay al menos un SArray, convertir todos a ese tipo.
    T1 = typeof(Σt)
    T2 = typeof(Σs0)
    if isnothing(χ)
        χ = zeros(T, g)
        χ[end] = one(T)
    else
        if !isone(sum(χ))
            error("χ no suma 1.")
        end
    end
    return CrossSections{g,T1,T2}(D, S, χ, Σt, Σa, νΣf, eΣf, Σs0, Σs1)
end