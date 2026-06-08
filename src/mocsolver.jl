import RayTracing:
    Track, Segment,
    DirectionType, Forward, Backward,
    universal_id, bc_fwd, bc_bwd, dir_next_track_fwd, dir_next_track_bwd

# TODO(feature): expose scalar flux as a Gridap field for interpolation and VTK output.
struct MoCSolution{T<:Real,P<:MoCProblem} <: TransportSolution
    prob::P

    keff::T
    residual::T
    iterations::Int

    # Region-major scalar flux, indexed by (region, group).
    φ::Vector{T}

    # Region-major reduced source, indexed by (region, group).
    q::Vector{T}

    # Region-integrated source used to compute the nonlinear iteration residual.
    Q::Vector{T}

    # Boundary angular flux for all tracks, directions, polar angles, and energy groups.
    boundary_ψ::Vector{T}
    start_boundary_ψ::Vector{T}
end

function MoCSolution{T}(prob::MoCProblem) where {T}
    @unpack nφ, nψ = prob
    NRegions = nregions(prob)
    φ = Vector{T}(undef, nφ)
    q = Vector{T}(undef, nφ)
    Q = Vector{T}(undef, NRegions)
    boundary_ψ = Vector{T}(undef, nψ)
    start_boundary_ψ = Vector{T}(undef, nψ)
    return MoCSolution(prob, one(T), zero(T), 0, φ, q, Q, boundary_ψ, start_boundary_ψ)
end

function show(io::IO, sol::MoCSolution)
    @unpack keff, residual, iterations = sol
    println(io, "  keff: ", keff)
    println(io, "  Residual: ", residual)
    print(io, "  Iterations: ", iterations)
end

function (sol::MoCSolution)(i::Int, g::Int)
    NGroups = ngroups(sol.prob)
    g in 1:NGroups || throw(DomainError(g, "`g` is outside of domain."))
    return sol.φ[@region_index(i, g)]
end

function (sol::MoCSolution)(g::Int)
    NGroups = ngroups(sol.prob)
    g in 1:NGroups || throw(DomainError(g, "`g` is outside of domain."))
    return view(sol.φ, g:NGroups:lastindex(sol.φ))
end

function solve(
    prob::MoCProblem;
    max_iterations::Int=1000, max_residual::Real=1e-7, debug::Bool=false
)
    return _solve_eigenvalue_problem(prob, max_iterations, max_residual, debug)
end

function _solve_eigenvalue_problem(prob::MoCProblem, max_iter::Int, max_ϵ::Real, debug::Bool)

    T = eltype(prob)
    sol = MoCSolution{T}(prob)

    optical_length!(prob)
    validate_track_volumes(prob)

    @set! sol.keff = one(T)
    set_uniform_φ!(sol, one(T))
    set_uniform_Q!(sol, zero(T))
    set_uniform_start_boundary_ψ!(sol, zero(T))
    update_boundary_ψ!(sol)

    debug && @info "MoC iterations start..."
    ϵ = Inf
    iter = 0
    while iter < max_iter

        normalize_fluxes!(sol, prob)
        compute_q!(sol, prob)
        compute_φ!(sol, prob)
        @set! sol.keff *= total_fission_source(sol, prob)
        ϵ = residual(sol, prob)

        debug && iszero(iter % 10) && @info "iteration $(iter)" sol.keff ϵ

        # The first two iterations seed boundary fluxes and the previous-source residual.
        # After that, `Q` contains a meaningful previous iterate.

        if iter > 1 && isless(ϵ, max_ϵ)
            break
        end

        iter += 1
    end

    @set! sol.residual = ϵ
    @set! sol.iterations = iter

    return sol
end

function optical_length!(prob::MoCProblem)
    @unpack trackgenerator = prob
    @unpack tracks_by_uid = trackgenerator

    for track in tracks_by_uid
        _optical_length!(prob, track)
    end

    return nothing
end

function _optical_length!(prob::MoCProblem, track::Track)
    NGroups = ngroups(prob)
    τ = optical_lengths!(prob, track)
    attenuation = attenuation_factors!(prob, track)
    @unpack sinθs = prob.quadrature.polar
    n_polar_half_count = n_polar_half(prob.quadrature.polar)

    for (s, segment) in enumerate(track.segments)
        @unpack ℓ = segment
        xs = getxs(prob, fsr_id(prob, segment.element))
        @unpack Σt = xs
        @inbounds for g in 1:NGroups
            τgs = Σt[g] * ℓ
            τ[g, s] = τgs
            for p in 1:n_polar_half_count
                attenuation[g, p, s] = -expm1(-τgs / sinθs[p])
            end
        end
    end

    return nothing
end

function optical_lengths!(prob::MoCProblem{Dim,NRegions,NGroups,T}, track::Track) where {Dim,NRegions,NGroups,T}
    @unpack optical_lengths = prob
    uid = universal_id(track)
    nsegments = length(track.segments)
    τ = optical_lengths[uid]

    if size(τ) != (NGroups, nsegments)
        τ = Matrix{T}(undef, NGroups, nsegments)
        optical_lengths[uid] = τ
    end

    return τ
end

function attenuation_factors!(
    prob::MoCProblem{Dim,NRegions,NGroups,T}, track::Track
) where {Dim,NRegions,NGroups,T}
    @unpack attenuation_factors, quadrature = prob
    uid = universal_id(track)
    nsegments = length(track.segments)
    n_polar_half_count = n_polar_half(quadrature.polar)
    factors = attenuation_factors[uid]

    if size(factors) != (NGroups, n_polar_half_count, nsegments)
        factors = Array{T,3}(undef, NGroups, n_polar_half_count, nsegments)
        attenuation_factors[uid] = factors
    end

    return factors
end

function validate_track_volumes(prob::MoCProblem)
    @unpack volumes = prob

    invalid = findall(v -> !isfinite(v) || v <= zero(v), volumes)
    if !isempty(invalid)
        n = length(invalid)
        sample = first(invalid, min(n, 5))
        throw(ArgumentError(
            "transport FSR volumes must be positive and finite; found $n " *
            "invalid volume(s), including region ids $(sample). Refine the track " *
            "spacing or enable volume correction before solving."
        ))
    end

    return nothing
end

@inline set_uniform_φ!(sol::MoCSolution, φ::Real) = fill!(sol.φ, φ)
@inline set_uniform_Q!(sol::MoCSolution, Q::Real) = fill!(sol.Q, Q)
@inline set_uniform_start_boundary_ψ!(sol::MoCSolution, ψ::Real) = fill!(sol.start_boundary_ψ, ψ)
@inline update_boundary_ψ!(sol::MoCSolution) = copy!(sol.boundary_ψ, sol.start_boundary_ψ)

function normalize_fluxes!(sol::MoCSolution, prob::MoCProblem)
    @unpack φ, boundary_ψ, start_boundary_ψ = sol

    qft = total_fission_source(sol, prob)

    λ = 1 / qft
    φ .*= λ
    boundary_ψ .*= λ
    start_boundary_ψ .*= λ

    return nothing
end

function total_fission_source(sol::MoCSolution{T}, prob::MoCProblem) where {T}
    NGroups = ngroups(prob)
    NRegions = nregions(prob)
    @unpack φ = sol
    @unpack volumes = prob

    qft = zero(T)
    @inbounds for i in 1:NRegions
        xs = getxs(prob, i)
        @unpack νΣf = xs
        fissionable = isfissionable(xs)

        if fissionable
            for g′ in 1:NGroups
                ig′ = @region_index(i, g′)
                qft += νΣf[g′] * φ[ig′] * volumes[i]
            end
        end
    end

    return qft
end

function compute_q!(sol::MoCSolution{T}, prob::MoCProblem) where {T}
    NGroups = ngroups(prob)
    NRegions = nregions(prob)
    @unpack keff, φ, q = sol

    @inbounds for i in 1:NRegions
        xs = getxs(prob, i)
        @unpack χ, Σt, νΣf, Σs0 = xs
        fissionable = isfissionable(xs)

        for g in 1:NGroups
            ig = @region_index(i, g)
            qig = zero(T)
            for g′ in 1:NGroups
                ig′ = @region_index(i, g′)
                qig += Σs0[g′, g] * φ[ig′]
                if fissionable
                    qig += 1 / keff * χ[g] * νΣf[g′] * φ[ig′]
                end
            end
            qig /= (4π * Σt[g])
            q[ig] = qig
        end
    end

    return nothing
end

function compute_φ!(sol::MoCSolution{T}, prob::MoCProblem) where {T}
    @unpack trackgenerator = prob
    @unpack tracks_by_uid = trackgenerator

    set_uniform_φ!(sol, zero(T))
    update_boundary_ψ!(sol)

    for track in tracks_by_uid
        tally!(sol, prob, track, Forward)
        tally!(sol, prob, track, Backward)
    end

    add_q_to_φ!(sol, prob)

    return nothing
end

function tally!(sol::MoCSolution, prob::MoCProblem, track::Track, dir::DirectionType)
    NGroups = ngroups(prob)
    @unpack quadrature = prob

    n_polar_half_count = n_polar_half(quadrature.polar)

    t = universal_id(track)
    d = Int32(dir)

    i = @angular_index(t, d, 1, 1)
    j = i + NGroups * n_polar_half_count - 1
    # Incoming boundary flux for one track and direction, indexed by (polar, group).
    boundary_ψ = @view sol.boundary_ψ[i:j]

    segments = track.segments
    attenuation = prob.attenuation_factors[t]

    if dir == Forward
        for s in eachindex(segments)
            tally_φ!(sol, prob, track, segments[s], attenuation, s, boundary_ψ)
        end
    elseif dir == Backward
        for s in lastindex(segments):-1:firstindex(segments)
            tally_φ!(sol, prob, track, segments[s], attenuation, s, boundary_ψ)
        end
    end

    set_start_boundary_ψ!(sol, prob, track, boundary_ψ, dir)

    return nothing
end

function tally_φ!(
    sol::MoCSolution,
    prob::MoCProblem,
    track::Track,
    segment::Segment,
    attenuation::AbstractArray{<:Real,3},
    segment_idx::Integer,
    boundary_ψ::AbstractVector
)
    NGroups = ngroups(prob)
    @unpack φ, q = sol
    @unpack quadrature = prob
    @unpack polar, ω = quadrature

    i = fsr_id(prob, segment.element)
    a = track.azim_idx
    n_polar_half_count = n_polar_half(polar)

    # TODO(performance): benchmark loop order for small NGroups versus larger polar sets.
    @inbounds for g in 1:NGroups, p in 1:n_polar_half_count
        pg = @reduced_angular_index(p, g)
        ig = @region_index(i, g)
        Δψ = (boundary_ψ[pg] - q[ig]) * attenuation[g, p, segment_idx]
        φ[ig] += 2 * ω[a, p] * Δψ  # Symmetry accounts for the omitted polar half-space.
        boundary_ψ[pg] -= Δψ
    end

    return nothing
end

function set_start_boundary_ψ!(
    sol::MoCSolution{T},
    prob::MoCProblem,
    current_track::Track,
    boundary_ψ::AbstractVector,
    dir::DirectionType
) where {T}

    NGroups = ngroups(prob)
    @unpack start_boundary_ψ = sol
    @unpack quadrature = prob

    n_polar_half_count = n_polar_half(quadrature.polar)

    if dir == Forward
        next_track = current_track.next_track_fwd
        next_track_dir = dir_next_track_fwd(current_track)
        flag = !isequal(bc_fwd(current_track), Vacuum)
    elseif dir == Backward
        next_track = current_track.next_track_bwd
        next_track_dir = dir_next_track_bwd(current_track)
        flag = !isequal(bc_bwd(current_track), Vacuum)
    end

    t = universal_id(next_track)
    d = Int32(next_track_dir)

    # TODO(performance): benchmark loop order against the boundary-flux memory layout.
    @inbounds for g in 1:NGroups, p in 1:n_polar_half_count
        tdpg = @angular_index(t, d, p, g)
        pg = @reduced_angular_index(p, g)
        start_boundary_ψ[tdpg] = flag ? boundary_ψ[pg] : zero(T)
    end

    return nothing
end

function add_q_to_φ!(sol::MoCSolution, prob::MoCProblem)
    NGroups = ngroups(prob)
    NRegions = nregions(prob)
    @unpack φ, q = sol
    @unpack volumes = prob

    @inbounds for i in 1:NRegions
        xs = getxs(prob, i)
        @unpack Σt = xs

        for g in 1:NGroups
            ig = @region_index(i, g)
            φ[ig] /= (Σt[g] * volumes[i])
            φ[ig] += (4π * q[ig])
        end
    end

    return nothing
end

function residual(sol::MoCSolution{T}, prob::MoCProblem) where {T}
    NGroups = ngroups(prob)
    NRegions = nregions(prob)
    @unpack keff, φ, Q = sol

    ϵ = zero(T)
    @inbounds for i in 1:NRegions

        xs = getxs(prob, i)
        @unpack νΣf, Σs0 = xs
        fissionable = isfissionable(xs)

        old_qi = Q[i]
        new_qi = zero(T)

        # Fission source contribution in this region. χ is normalized, so it does not enter
        # the group-integrated source used for residual convergence.
        if fissionable
            for g′ in 1:NGroups
                ig′ = @region_index(i, g′)
                νΣfg′ = νΣf[g′]
                new_qi += νΣfg′ * φ[ig′]
            end
        end

        new_qi /= keff

        # Scattering source contribution in this region.
        for g in 1:NGroups, g′ in 1:NGroups
            ig′ = @region_index(i, g′)
            Σsgg′ = Σs0[g′, g]
            new_qi += Σsgg′ * φ[ig′]
        end

        if old_qi > 0
            ϵ += ((new_qi - old_qi) / old_qi)^2
        end

        Q[i] = new_qi
    end

    ϵ = sqrt(ϵ / NRegions)

    return ϵ
end
