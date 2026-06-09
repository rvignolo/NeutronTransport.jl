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

    # Reused by opt-in threaded kernels to avoid per-iteration scratch allocation.
    thread_sums::Vector{T}
    thread_φ::Vector{Vector{T}}
    parallel_sweep_checked::Vector{Bool}
end

function MoCSolution{T}(prob::MoCProblem) where {T}
    @unpack nφ, nψ = prob
    NRegions = nregions(prob)
    φ = Vector{T}(undef, nφ)
    q = Vector{T}(undef, nφ)
    Q = Vector{T}(undef, NRegions)
    boundary_ψ = Vector{T}(undef, nψ)
    start_boundary_ψ = Vector{T}(undef, nψ)
    thread_sums = Vector{T}(undef, Base.Threads.nthreads())
    thread_φ = Vector{Vector{T}}()
    parallel_sweep_checked = [false]
    return MoCSolution(
        prob, one(T), zero(T), 0, φ, q, Q, boundary_ψ, start_boundary_ψ,
        thread_sums, thread_φ, parallel_sweep_checked
    )
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
    max_iterations::Int=1000,
    max_residual::Real=1e-7,
    debug::Bool=false,
    parallel::Bool=false
)
    return _solve_eigenvalue_problem(prob, max_iterations, max_residual, debug, parallel)
end

function _solve_eigenvalue_problem(
    prob::MoCProblem, max_iter::Int, max_ϵ::Real, debug::Bool, parallel::Bool
)

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
    max_ϵ = T(max_ϵ)
    ϵ = T(Inf)
    iter = 0
    while iter < max_iter

        normalize_fluxes!(sol, prob; parallel=parallel)
        compute_q!(sol, prob; parallel=parallel)
        compute_φ!(sol, prob; parallel=parallel)
        @set! sol.keff *= total_fission_source(sol, prob; parallel=parallel)
        ϵ = residual(sol, prob; parallel=parallel)

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
@inline use_threads(parallel::Bool) = parallel && Base.Threads.nthreads() > 1

function thread_sums!(sol::MoCSolution{T}) where {T}
    n = Base.Threads.nthreads()
    partials = sol.thread_sums
    if length(partials) != n
        resize!(partials, n)
    end
    fill!(partials, zero(T))
    return partials
end

function thread_φ_buffers!(sol::MoCSolution{T}) where {T}
    n = Base.Threads.nthreads()
    nφ = length(sol.φ)
    buffers = sol.thread_φ

    if length(buffers) != n
        resize!(buffers, n)
    end

    for tid in 1:n
        if !isassigned(buffers, tid) || length(buffers[tid]) != nφ
            buffers[tid] = zeros(T, nφ)
        else
            fill!(buffers[tid], zero(T))
        end
    end

    return buffers
end

function normalize_fluxes!(sol::MoCSolution, prob::MoCProblem; parallel::Bool=false)
    @unpack φ, boundary_ψ, start_boundary_ψ = sol

    qft = total_fission_source(sol, prob; parallel=parallel)

    λ = one(qft) / qft
    φ .*= λ
    boundary_ψ .*= λ
    start_boundary_ψ .*= λ

    return nothing
end

function total_fission_source(
    sol::MoCSolution{T}, prob::MoCProblem; parallel::Bool=false
) where {T}
    if use_threads(parallel)
        return total_fission_source_threaded(sol, prob)
    end

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

function total_fission_source_threaded(sol::MoCSolution{T}, prob::MoCProblem) where {T}
    NGroups = ngroups(prob)
    NRegions = nregions(prob)
    nworkers = Base.Threads.nthreads()
    @unpack φ = sol
    @unpack volumes = prob

    partials = thread_sums!(sol)

    Base.Threads.@threads for part in 1:nworkers
        qft_part = zero(T)
        @inbounds for i in part:nworkers:NRegions
            xs = getxs(prob, i)
            if isfissionable(xs)
                @unpack νΣf = xs
                volume = volumes[i]
                for g′ in 1:NGroups
                    ig′ = @region_index(i, g′)
                    qft_part += νΣf[g′] * φ[ig′] * volume
                end
            end
        end
        partials[part] = qft_part
    end

    qft = zero(T)
    @inbounds for part in 1:nworkers
        qft += partials[part]
    end

    return qft
end

function compute_q!(sol::MoCSolution{T}, prob::MoCProblem; parallel::Bool=false) where {T}
    NRegions = nregions(prob)
    @unpack keff, φ, q = sol
    inv_keff = inv(keff)

    if use_threads(parallel)
        Base.Threads.@threads for i in 1:NRegions
            compute_q_region!(q, φ, prob, i, inv_keff)
        end
        return nothing
    end

    @inbounds for i in 1:NRegions
        compute_q_region!(q, φ, prob, i, inv_keff)
    end

    return nothing
end

function compute_q_region!(
    q::Vector{T},
    φ::Vector{T},
    prob::MoCProblem{Dim,NRegions,NGroups,T},
    i::Int,
    inv_keff::T
) where {T,Dim,NRegions,NGroups}
    xs = getxs(prob, i)
    @unpack χ, Σt, νΣf, Σs0 = xs
    fissionable = isfissionable(xs)
    inv_4π = inv(T(4π))

    fission_source = zero(T)
    if fissionable
        @inbounds for g′ in 1:NGroups
            ig′ = @region_index(i, g′)
            fission_source += νΣf[g′] * φ[ig′]
        end
    end

    @inbounds for g in 1:NGroups
        ig = @region_index(i, g)
        scattering_source = zero(T)
        for g′ in 1:NGroups
            ig′ = @region_index(i, g′)
            scattering_source += Σs0[g′, g] * φ[ig′]
        end
        fission_term = fissionable ? inv_keff * χ[g] * fission_source : zero(T)
        q[ig] = (scattering_source + fission_term) * inv_4π / Σt[g]
    end

    return nothing
end

function compute_φ!(
    sol::MoCSolution{T}, prob::MoCProblem; parallel::Bool=false
) where {T}
    @unpack trackgenerator = prob
    @unpack tracks_by_uid = trackgenerator

    set_uniform_φ!(sol, zero(T))
    update_boundary_ψ!(sol)

    if use_threads(parallel)
        compute_φ_threaded!(sol, prob, tracks_by_uid)
    else
        for track in tracks_by_uid
            tally!(sol, prob, track, Forward)
            tally!(sol, prob, track, Backward)
        end
        add_q_to_φ!(sol, prob; parallel=false)
    end

    return nothing
end

function compute_φ_threaded!(
    sol::MoCSolution{T}, prob::MoCProblem, tracks_by_uid
) where {T}
    ensure_parallel_sweep_supported!(sol, prob)
    φ_buffers = thread_φ_buffers!(sol)
    fill!(sol.start_boundary_ψ, zero(T))

    nworkers = length(φ_buffers)
    Base.@sync for worker in 1:nworkers
        Base.Threads.@spawn begin
            φ_thread = φ_buffers[worker]
            for track_idx in worker:nworkers:length(tracks_by_uid)
                track = tracks_by_uid[track_idx]
                tally!(sol, prob, track, Forward, φ_thread; write_vacuum=false)
                tally!(sol, prob, track, Backward, φ_thread; write_vacuum=false)
            end
        end
    end

    reduce_thread_φ_and_add_q!(sol, prob, φ_buffers)

    return nothing
end

function reduce_thread_φ_and_add_q!(
    sol::MoCSolution{T},
    prob::MoCProblem{Dim,NRegions,NGroups,T},
    φ_buffers::Vector{Vector{T}},
) where {T,Dim,NRegions,NGroups}
    @unpack φ, q = sol
    @unpack volumes = prob
    nworkers = length(φ_buffers)
    fourπ = T(4π)

    Base.Threads.@threads for i in 1:NRegions
        xs = getxs(prob, i)
        @unpack Σt = xs
        volume = volumes[i]
        @inbounds for g in 1:NGroups
            ig = @region_index(i, g)
            φig = zero(T)
            for tid in 1:nworkers
                φig += φ_buffers[tid][ig]
            end
            φ[ig] = φig / (Σt[g] * volume) + fourπ * q[ig]
        end
    end

    return nothing
end

function ensure_parallel_sweep_supported!(sol::MoCSolution, prob::MoCProblem)
    if !sol.parallel_sweep_checked[1]
        validate_parallel_sweep_boundaries(prob)
        sol.parallel_sweep_checked[1] = true
    end
    return nothing
end

function validate_parallel_sweep_boundaries(prob::MoCProblem)
    @unpack tracks_by_uid = prob.trackgenerator
    target_written = falses(2 * length(tracks_by_uid))

    for track in tracks_by_uid
        validate_parallel_boundary_writer!(target_written, track, Forward)
        validate_parallel_boundary_writer!(target_written, track, Backward)
    end

    return nothing
end

function validate_parallel_boundary_writer!(
    target_written::AbstractVector{Bool}, current_track::Track, dir::DirectionType
)
    if dir == Forward
        isequal(bc_fwd(current_track), Vacuum) && return nothing
        next_track = current_track.next_track_fwd
        next_track_dir = dir_next_track_fwd(current_track)
    elseif dir == Backward
        isequal(bc_bwd(current_track), Vacuum) && return nothing
        next_track = current_track.next_track_bwd
        next_track_dir = dir_next_track_bwd(current_track)
    else
        throw(ArgumentError("unsupported sweep direction `$dir`."))
    end

    target_track_id = universal_id(next_track)
    target_dir = Int(next_track_dir)
    target_dir in (0, 1) ||
        throw(ArgumentError("parallel sweep expected zero-based track directions."))
    target_idx = 2 * (target_track_id - 1) + target_dir + 1
    target_idx in eachindex(target_written) ||
        throw(ArgumentError("parallel sweep boundary target is outside track storage."))

    if target_written[target_idx]
        throw(ArgumentError(
            "parallel sweep requires unique non-vacuum boundary writers; " *
            "multiple outgoing tracks target track $target_track_id, direction $target_dir."
        ))
    end

    target_written[target_idx] = true

    return nothing
end

function tally!(sol::MoCSolution, prob::MoCProblem, track::Track, dir::DirectionType)
    return tally!(sol, prob, track, dir, sol.φ; write_vacuum=true)
end

function tally!(
    sol::MoCSolution{T},
    prob::MoCProblem{Dim,NRegions,NGroups,T},
    track::Track{T},
    dir::DirectionType,
    φ_accum::Vector{T};
    write_vacuum::Bool=true
) where {T,Dim,NRegions,NGroups}
    @unpack quadrature = prob

    n_polar_half_count = n_polar_half(quadrature.polar)

    t = universal_id(track)
    d = Int(dir)

    # Offset into incoming boundary flux for one track and direction, indexed by
    # contiguous (polar, group) entries. Keeping an integer offset avoids creating a
    # SubArray for every track/direction sweep.
    boundary_offset = @angular_index(t, d, 1, 1) - 1
    boundary_ψ = sol.boundary_ψ

    segments = track.segments
    attenuation = prob.attenuation_factors[t]

    if dir == Forward
        for s in eachindex(segments)
            tally_φ!(
                sol, prob, track, segments[s], attenuation, s,
                boundary_ψ, boundary_offset, φ_accum
            )
        end
    elseif dir == Backward
        for s in lastindex(segments):-1:firstindex(segments)
            tally_φ!(
                sol, prob, track, segments[s], attenuation, s,
                boundary_ψ, boundary_offset, φ_accum
            )
        end
    else
        throw(ArgumentError("unsupported sweep direction `$dir`."))
    end

    set_start_boundary_ψ!(
        sol, prob, track, boundary_ψ, boundary_offset, dir; write_vacuum=write_vacuum
    )

    return nothing
end

function tally_φ!(
    sol::MoCSolution{T},
    prob::MoCProblem{Dim,NRegions,NGroups,T},
    track::Track{T},
    segment::Segment{T},
    attenuation::Array{T,3},
    segment_idx::Int,
    boundary_ψ::Vector{T},
    boundary_offset::Int,
    φ_accum::Vector{T}
) where {T,Dim,NRegions,NGroups}
    @unpack q = sol
    @unpack quadrature = prob
    @unpack polar, ω = quadrature

    i = Int(fsr_id(prob, segment.element))
    a = track.azim_idx
    n_polar_half_count = n_polar_half(polar)
    region_offset = (i - 1) * NGroups

    @inbounds for p in 1:n_polar_half_count
        polar_offset = (p - 1) * NGroups
        weight = 2 * ω[a, p]  # Symmetry accounts for omitted polar half-space.
        for g in 1:NGroups
            ψ_idx = boundary_offset + polar_offset + g
            ig = region_offset + g
            Δψ = (boundary_ψ[ψ_idx] - q[ig]) * attenuation[g, p, segment_idx]
            φ_accum[ig] += weight * Δψ
            boundary_ψ[ψ_idx] -= Δψ
        end
    end

    return nothing
end

function set_start_boundary_ψ!(
    sol::MoCSolution{T},
    prob::MoCProblem{Dim,NRegions,NGroups,T},
    current_track::Track{T},
    boundary_ψ::Vector{T},
    boundary_offset::Int,
    dir::DirectionType;
    write_vacuum::Bool=true
) where {T,Dim,NRegions,NGroups}
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
    else
        throw(ArgumentError("unsupported sweep direction `$dir`."))
    end

    if !flag && !write_vacuum
        return nothing
    end

    t = universal_id(next_track)
    d = Int32(next_track_dir)

    @inbounds for p in 1:n_polar_half_count
        polar_offset = (p - 1) * NGroups
        for g in 1:NGroups
            tdpg = @angular_index(t, d, p, g)
            ψ_idx = boundary_offset + polar_offset + g
            start_boundary_ψ[tdpg] = flag ? boundary_ψ[ψ_idx] : zero(T)
        end
    end

    return nothing
end

function add_q_to_φ!(sol::MoCSolution, prob::MoCProblem; parallel::Bool=false)
    NRegions = nregions(prob)
    @unpack φ, q = sol
    @unpack volumes = prob

    if use_threads(parallel)
        Base.Threads.@threads for i in 1:NRegions
            add_q_region_to_φ!(φ, q, prob, volumes, i)
        end
        return nothing
    end

    @inbounds for i in 1:NRegions
        add_q_region_to_φ!(φ, q, prob, volumes, i)
    end

    return nothing
end

function add_q_region_to_φ!(
    φ::Vector{T},
    q::Vector{T},
    prob::MoCProblem{Dim,NRegions,NGroups,T},
    volumes::Vector{T},
    i::Int
) where {T,Dim,NRegions,NGroups}
    xs = getxs(prob, i)
    @unpack Σt = xs
    fourπ = T(4π)

    @inbounds for g in 1:NGroups
        ig = @region_index(i, g)
        φ[ig] /= (Σt[g] * volumes[i])
        φ[ig] += fourπ * q[ig]
    end

    return nothing
end

function residual(sol::MoCSolution{T}, prob::MoCProblem; parallel::Bool=false) where {T}
    if use_threads(parallel)
        return residual_threaded(sol, prob)
    end

    NRegions = nregions(prob)
    @unpack Q = sol

    ϵ = zero(T)
    @inbounds for i in 1:NRegions
        old_qi = Q[i]
        new_qi = region_source(sol, prob, i)

        if old_qi > 0
            ϵ += ((new_qi - old_qi) / old_qi)^2
        end

        Q[i] = new_qi
    end

    ϵ = sqrt(ϵ / NRegions)

    return ϵ
end

function residual_threaded(sol::MoCSolution{T}, prob::MoCProblem) where {T}
    NRegions = nregions(prob)
    nworkers = Base.Threads.nthreads()
    @unpack Q = sol

    partials = thread_sums!(sol)

    Base.Threads.@threads for part in 1:nworkers
        ϵ_part = zero(T)
        @inbounds for i in part:nworkers:NRegions
            new_qi = region_source(sol, prob, i)
            old_qi = Q[i]

            if old_qi > 0
                ϵ_part += ((new_qi - old_qi) / old_qi)^2
            end

            Q[i] = new_qi
        end
        partials[part] = ϵ_part
    end

    ϵ = zero(T)
    @inbounds for part in 1:nworkers
        ϵ += partials[part]
    end

    return sqrt(ϵ / NRegions)
end

function region_source(
    sol::MoCSolution{T}, prob::MoCProblem{Dim,NRegions,NGroups}, i::Integer
) where {T,Dim,NRegions,NGroups}
    @unpack keff, φ = sol
    xs = getxs(prob, i)
    @unpack νΣf, Σs0_sum = xs

    new_qi = zero(T)

    # Fission source contribution in this region. χ is normalized, so it does not enter
    # the group-integrated source used for residual convergence.
    if isfissionable(xs)
        @inbounds for g′ in 1:NGroups
            ig′ = @region_index(i, g′)
            new_qi += νΣf[g′] * φ[ig′]
        end
    end

    new_qi /= keff

    # Group-integrated scattering source contribution in this region.
    @inbounds for g′ in 1:NGroups
        ig′ = @region_index(i, g′)
        new_qi += Σs0_sum[g′] * φ[ig′]
    end

    return new_qi
end
