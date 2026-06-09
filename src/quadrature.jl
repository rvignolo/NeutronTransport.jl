import RayTracing: AzimuthalQuadrature, n_azim_half

"""
    Quadrature{A<:AzimuthalQuadrature,P<:PolarQuadrature,T<:Real}

Holds information regarding both the azimuthal and polar quadrature and the total weights.
"""
struct Quadrature{T<:Real,A<:AzimuthalQuadrature,P<:PolarQuadrature}
    azimuthal::A
    polar::P
    ω::Matrix{T}
end

function Quadrature(
    azimuthal::AzimuthalQuadrature{Na,N2,N4,T}, polar::PolarQuadrature{Np,T}
) where {Na,N2,N4,Np,T}
    @unpack δs, ωₐ = azimuthal
    @unpack sinθs, ωₚ = polar
    n_azim_half_count = n_azim_half(azimuthal)
    n_polar_half_count = n_polar_half(polar)

    # TODO(performance): store only quadrant weights if the sweep can reuse symmetry.
    ω = Matrix{T}(undef, n_azim_half_count, n_polar_half_count)

    for i in 1:n_azim_half_count, j in 1:n_polar_half_count
        ω[i, j] = 4π * ωₐ[i] * ωₚ[j] * δs[i] * sinθs[j]
    end
    return Quadrature(azimuthal, polar, ω)
end
