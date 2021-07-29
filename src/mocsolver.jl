import RayTracing:
    Track, Segment,
    DirectionType, Forward, Backward,
    universal_id, bc_fwd, bc_bwd, dir_next_track_fwd, dir_next_track_bwd

# TODO: Gridap solution object might be useful for interpolation at any point in space...
struct MoCSolution{T<:Real,P<:MoCProblem} <: TransportSolution
    prob::P

    keff::T
    residual::T
    iterations::Int

    # scalar flux
    φ::Vector{T}

    # reduced source
    q::Vector{T}

    # total source (integrated in energy given per region) used for convergence purposes
    Q::Vector{T}

    # angular flux at boundary
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
    print(io,   "  Iterations: ", iterations)
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

    T =  eltype(prob)
    sol = MoCSolution{T}(prob)

    optical_length!(prob)

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

        # we need at least three iterations:
        #  0. boundary conditions do not exists (unless everything is vaccum)
        #  1. we can compute a residual but it is not correct (previous iteration is fruit)
        #  2. now we can compute a correct residual

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

    for segment in track.segments
        @unpack ℓ, τ = segment

        resize!(τ, NGroups)

        xs = getxs(prob, segment.element)
        @unpack Σt = xs
        @inbounds for g in 1:NGroups
            τ[g] = Σt[g] * ℓ
        end
    end
end

@inline set_uniform_φ!(sol::MoCSolution, φ::Real) = fill!(sol.φ, φ)
@inline set_uniform_Q!(sol::MoCSolution, Q::Real) = fill!(sol.Q, Q)
@inline set_uniform_start_boundary_ψ!(sol::MoCSolution, ψ::Real) = fill!(sol.start_boundary_ψ, ψ)
@inline update_boundary_ψ!(sol::MoCSolution) = copy!(sol.boundary_ψ, sol.start_boundary_ψ)

function normalize_fluxes!(sol::MoCSolution, prob::MoCProblem)
    @unpack φ, boundary_ψ, start_boundary_ψ = sol

    # total fission source
    qft = total_fission_source(sol, prob)

    # λ is the normalization factor
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
    @unpack trackgenerator = prob
    @unpack volumes = trackgenerator

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
            # qft += sum(
            #     νΣf[g′] * φ[@region_index(i, g′)] * volumes[i] for g′ in 1:NGroups
            # )
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
            # qig = sum(
            #     Σs0[g′, g] * φ[@region_index(i, g′)] +
            #     1 / keff * χ[g] * νΣf[g′] * φ[@region_index(i, g′)]
            #     for g′ in 1:NGroups
            # )

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
    update_boundary_ψ!(sol) # update entry ψ for all rays (including d, p and g dependence)

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

    n_polar_2 = npolar2(quadrature.polar)

    t = universal_id(track)
    d = Int32(dir)

    i = @angular_index(t, d, 1, 1)
    j = i + NGroups * n_polar_2 - 1
    boundary_ψ = @view sol.boundary_ψ[i:j] # boundary_ψ in for a given track in a given direction as function of (p, g)

    segments = dir == Forward ? track.segments : reverse!(track.segments)

    for segment in segments
        tally_φ!(sol, prob, track, segment, boundary_ψ)
    end

    set_start_boundary_ψ!(sol, prob, track, boundary_ψ, dir)

    if dir == Backward
        reverse!(track.segments)
    end

    return nothing
end

function tally_φ!(
    sol::MoCSolution,
    prob::MoCProblem,
    track::Track,
    segment::Segment,
    boundary_ψ::AbstractVector
)
    NGroups = ngroups(prob)
    @unpack φ, q = sol
    @unpack quadrature = prob
    @unpack polar, ω = quadrature
    @unpack sinθs = polar
    @unpack τ = segment

    i = segment.element
    a = track.azim_idx
    n_polar_2 = npolar2(polar)

    # TODO: is this the best possible loop order?
    @inbounds for g in 1:NGroups, p in 1:n_polar_2
        pg = @reduced_angular_index(p, g)
        ig = @region_index(i, g)
        Δψ = (boundary_ψ[pg] - q[ig]) * (1 - exp(-τ[g] / sinθs[p]))
        φ[ig] += 2 * ω[a, p] * Δψ # multiplied by 2 because we only loop in n_polar_2
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

    n_polar_2 = npolar2(quadrature.polar)

    if dir == Forward
        next_track = current_track.next_track_fwd
        next_track_dir = dir_next_track_fwd(current_track)
        flag = !isequal(bc_fwd(current_track), Vaccum)
    elseif dir == Backward
        next_track = current_track.next_track_bwd
        next_track_dir = dir_next_track_bwd(current_track)
        flag = !isequal(bc_bwd(current_track), Vaccum)
    end

    t = universal_id(next_track)
    d = Int32(next_track_dir)

    # TODO: is this the best possible loop order?
    @inbounds for g in 1:NGroups, p in 1:n_polar_2
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
    @unpack trackgenerator = prob
    @unpack volumes = trackgenerator

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

        # total fission source in each region (χ sum 1 when ∫ in g)
        if fissionable
            for g′ in 1:NGroups
                ig′ = @region_index(i, g′)
                νΣfg′ = νΣf[g′]
                new_qi += νΣfg′ * φ[ig′]
            end
            # new_qi = sum(νΣf[g′] * φ[ig′] for g′ in 1:NGroups)
        end

        new_qi /= keff

        # total scattering source
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
