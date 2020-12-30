"""
    Quadrature{A<:AzimuthalQuadrature,P<:PolarQuadrature,T<:Real}

Holds information regarding both the azimuthal and polar quadrature and the total weights.
"""
struct Quadrature{A<:AzimuthalQuadrature,P<:PolarQuadrature,T<:Real}
    azimuthal::A
    polar::P
    ω::Matrix{T}
end

function Quadrature(azimuthal::AzimuthalQuadrature{T}, polar::PolarQuadrature{N,T}) where {N,T}
    @unpack n_azim_2, δs, ωₐ = azimuthal
    @unpack sinθs, ωₚ = polar
    n_polar_2 = npolar2(polar)

    # IDEA: I think we can use n_azim_4 because the matrix has repeated values
    ω = Matrix{T}(undef, n_azim_2, n_polar_2) # we have polar symmetry

    for i in RayTracing.both_dir(azimuthal), j in 1:n_polar_2
        ω[i, j] = 4π * ωₐ[i] * ωₚ[j] * δs[i] * sinθs[j]
    end
    # ω .*= 4π

    return Quadrature(azimuthal, polar, ω)
end