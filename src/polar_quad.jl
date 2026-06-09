
"""
    PolarQuadrature{N,T<:Real}

Holds information of the polar quadrature, such as the polar angles `θs` and the polar
weights `ωₚ`.
"""
struct PolarQuadrature{N,T<:Real}
    sinθs::Vector{T}  # TODO(performance): store quadrature data as SVectors.
    θs::Vector{T}
    ωₚ::Vector{T}
end

n_polar(::PolarQuadrature{N}) where {N} = N
n_polar_half(::PolarQuadrature{N}) where {N} = div(N, 2)

for polar_quad in (:TabuchiYamamoto, :GaussLegendre, :EqualWeight, :EqualAngle, :Leonard)
    @eval begin
        $polar_quad(N::Int, T::Type{<:Real}=Float64) = $polar_quad(Val(N), T)
    end
end

const _TL_NΘ_ = (2, 4, 6)
const _GL_NΘ_ = (2, 3, 4, 6, 8, 10, 12)

function TabuchiYamamoto(n_polar::Val{N}, T::Type{<:Real}=Float64) where {N}
    N in _TL_NΘ_ || throw(DomainError(N, "number of polar angles for TabuchiYamamoto must" *
                                         " be in $(string(_TL_NΘ_))."))

    sinθs, θs, ωₚ = ntuple(_ -> Vector{T}(undef, N), 3)

    # Nuclear Engineering Handbook, p. 1154.
    set_tabuchiyamamoto_data!(n_polar, sinθs, ωₚ)

    # Mirror the tabulated half-range values into the supplementary polar angles.
    for i in 1:div(N, 2)
        θs[i] = asin(sinθs[i])

        j = N - i + 1
        sinθs[j] = sinθs[i]
        θs[j] = π - θs[i]
        ωₚ[j] = ωₚ[i]
    end

    return PolarQuadrature{N,T}(sinθs, θs, ωₚ)
end

function set_tabuchiyamamoto_data!(::Val{2}, sinθs, ωₚ)
    sinθs[1] = 0.798184
    ωₚ[1] = 1 / 2
    return nothing
end

function set_tabuchiyamamoto_data!(::Val{4}, sinθs, ωₚ)
    sinθs[1] = 0.363900
    sinθs[2] = 0.899900
    ωₚ[1] = 0.212854 / 2
    ωₚ[2] = 0.787146 / 2
    return nothing
end

function set_tabuchiyamamoto_data!(::Val{6}, sinθs, ωₚ)
    sinθs[1] = 0.166648
    sinθs[2] = 0.537707
    sinθs[3] = 0.932954
    ωₚ[1] = 0.046233 / 2
    ωₚ[2] = 0.283619 / 2
    ωₚ[3] = 0.670148 / 2
    return nothing
end

function GaussLegendre(n_polar::Val{N}, T::Type{<:Real}=Float64) where {N}
    N in _GL_NΘ_ || throw(DomainError(N, "number of polar angles for GaussLegendre must be in $(string(_GL_NΘ_))."))

    sinθs, θs, ωₚ = ntuple(_ -> Vector{T}(undef, N), 3)

    # Nuclear Engineering Handbook, p. 1153.
    set_gausslegendre_data!(n_polar, θs, ωₚ)

    # Mirror the tabulated half-range values into the supplementary polar angles.
    for i in 1:div(N, 2)
        sinθs[i] = sin(θs[i])

        j = N - i + 1
        sinθs[j] = sinθs[i]
        θs[j] = π - θs[i]
        ωₚ[j] = ωₚ[i]
    end

    return PolarQuadrature{N,T}(sinθs, θs, ωₚ)
end

function set_gausslegendre_data!(::Val{2}, θs, ωₚ)
    θs[1] = acos(0.5773502691)
    ωₚ[1] = 1 / 2
    return nothing
end

function set_gausslegendre_data!(::Val{4}, θs, ωₚ)
    θs[1] = acos(0.3399810435)
    θs[2] = acos(0.8611363115)
    ωₚ[1] = 0.6521451549 / 2
    ωₚ[2] = 0.3478548451 / 2
    return nothing
end

function set_gausslegendre_data!(::Val{6}, θs, ωₚ)
    θs[1] = acos(0.2386191860)
    θs[2] = acos(0.6612093864)
    θs[3] = acos(0.9324695142)
    ωₚ[1] = 0.4679139346 / 2
    ωₚ[2] = 0.3607615730 / 2
    ωₚ[3] = 0.1713244924 / 2
    return nothing
end

function set_gausslegendre_data!(::Val{8}, θs, ωₚ)
    θs[1] = acos(0.1834346424)
    θs[2] = acos(0.5255324099)
    θs[3] = acos(0.7966664774)
    θs[4] = acos(0.9602898564)
    ωₚ[1] = 0.3626837834 / 2
    ωₚ[2] = 0.3137066459 / 2
    ωₚ[3] = 0.2223810344 / 2
    ωₚ[4] = 0.1012285363 / 2
    return nothing
end

function set_gausslegendre_data!(::Val{10}, θs, ωₚ)
    θs[1] = acos(0.1488743387)
    θs[2] = acos(0.4333953941)
    θs[3] = acos(0.6794095682)
    θs[4] = acos(0.8650633666)
    θs[5] = acos(0.9739065285)
    ωₚ[1] = 0.2955242247 / 2
    ωₚ[2] = 0.2692667193 / 2
    ωₚ[3] = 0.2190863625 / 2
    ωₚ[4] = 0.1494513492 / 2
    ωₚ[5] = 0.0666713443 / 2
    return nothing
end

function set_gausslegendre_data!(::Val{12}, θs, ωₚ)
    θs[1] = acos(0.1252334085)
    θs[2] = acos(0.3678314989)
    θs[3] = acos(0.5873179542)
    θs[4] = acos(0.7699026741)
    θs[5] = acos(0.9041172563)
    θs[6] = acos(0.9815606342)
    ωₚ[1] = 0.2491470458 / 2
    ωₚ[2] = 0.2334925365 / 2
    ωₚ[3] = 0.2031674267 / 2
    ωₚ[4] = 0.1600783286 / 2
    ωₚ[5] = 0.1069393260 / 2
    ωₚ[6] = 0.0471753364 / 2
    return nothing
end

function EqualWeight(n_polar::Val{N}, T::Type{<:Real}=Float64) where {N}
    iseven(N) || throw(DomainError(N, "number of polar angles for EqualWeight must be even."))

    sinθs, θs, ωₚ = ntuple(_ -> Vector{T}(undef, N), 3)
    half = div(N, 2)
    cosθa = one(T)

    for i in 1:half
        cosθb = cosθa - inv(T(half))
        θs[i] = acos((cosθa + cosθb) / 2)
        sinθs[i] = sin(θs[i])
        cosθa = cosθb
    end

    set_equal_spacing_weights!(ωₚ, θs, Val(N))
    mirror_polar_half!(sinθs, θs, ωₚ, Val(N))

    return PolarQuadrature{N,T}(sinθs, θs, ωₚ)
end

function EqualAngle(n_polar::Val{N}, T::Type{<:Real}=Float64) where {N}
    iseven(N) || throw(DomainError(N, "number of polar angles for EqualAngle must be even."))

    sinθs, θs, ωₚ = ntuple(_ -> Vector{T}(undef, N), 3)
    half = div(N, 2)
    Δθ = T(π) / T(N)
    θa = zero(T)

    for i in 1:half
        θb = θa + Δθ
        θs[i] = acos((cos(θa) + cos(θb)) / 2)
        sinθs[i] = sin(θs[i])
        θa = θb
    end

    set_equal_spacing_weights!(ωₚ, θs, Val(N))
    mirror_polar_half!(sinθs, θs, ωₚ, Val(N))

    return PolarQuadrature{N,T}(sinθs, θs, ωₚ)
end

function set_equal_spacing_weights!(ωₚ, θs, ::Val{N}) where {N}
    half = div(N, 2)

    for i in 1:half
        y1 = if i < half
            (cos(θs[i]) - cos(θs[i+1])) / 2
        else
            cos(θs[i])
        end

        y2 = if i > 1
            (cos(θs[i-1]) - cos(θs[i])) / 2
        else
            1 - cos(θs[i])
        end

        ωₚ[i] = (y1 + y2) / 2
    end

    return nothing
end

function mirror_polar_half!(sinθs, θs, ωₚ, ::Val{N}) where {N}
    for i in 1:div(N, 2)
        j = N - i + 1
        sinθs[j] = sinθs[i]
        θs[j] = π - θs[i]
        ωₚ[j] = ωₚ[i]
    end

    return nothing
end

function Leonard(n_polar::Val{N}, T::Type{<:Real}=Float64) where {N}
    N in _TL_NΘ_ || throw(DomainError(N, "number of polar angles for Leonard must be in $(string(_TL_NΘ_))."))

    sinθs, θs, ωₚ = ntuple(_ -> Vector{T}(undef, N), 3)

    # TODO(numerics): choose one authoritative Leonard quadrature source. Tabuchi,
    # OpenMOC, and the Hebert values in the Nuclear Engineering Handbook differ.
    set_leonard_data!(n_polar, sinθs, ωₚ)

    # Mirror the tabulated half-range values into the supplementary polar angles.
    for i in 1:div(N, 2)
        θs[i] = asin(sinθs[i])

        j = N - i + 1
        sinθs[j] = sinθs[i]
        θs[j] = π - θs[i]
        ωₚ[j] = ωₚ[i]
    end

    return PolarQuadrature{N,T}(sinθs, θs, ωₚ)
end

function set_leonard_data!(::Val{2}, sinθs, ωₚ)
    sinθs[1] = 0.752244
    ωₚ[1] = 1 / 2
    return nothing
end

function set_leonard_data!(::Val{4}, sinθs, ωₚ)
    sinθs[1] = 0.273658
    sinθs[2] = 0.865714
    ωₚ[1] = 0.139473 / 2
    ωₚ[2] = 0.860527 / 2
    return nothing
end

function set_leonard_data!(::Val{6}, sinθs, ωₚ)
    sinθs[1] = 0.103840
    sinθs[2] = 0.430723
    sinθs[3] = 0.905435
    ωₚ[1] = 0.020530 / 2
    ωₚ[2] = 0.219161 / 2
    ωₚ[3] = 0.760309 / 2
    return nothing
end
