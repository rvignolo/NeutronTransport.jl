import RayTracing:
    Track, Segment,
    DirectionType, Forward, Backward,
    universal_id, dir_next_track_fwd, dir_next_track_bwd, boundary_in, boundary_out

# TODO: Gridap solution object might be useful for interpolation at any point in space...
struct MoCSolution{T<:Real} <: TransportSolution
    keff::T
    residual::T
    iterations::Int

    # scalar flux
    φ::Vector{T}
    φ_prev::Vector{T}

    # reduced source
    q::Vector{T}

    # angular flux at boundary
    boundary_ψ::Vector{T}
    start_boundary_ψ::Vector{T}
end

function MoCSolution{T}(nφ::Int, nψ::Int) where {T}
    φ = Vector{T}(undef, nφ)
    φ_prev = Vector{T}(undef, nφ)
    q = Vector{T}(undef, nφ)
    boundary_ψ = Vector{T}(undef, nψ)
    start_boundary_ψ = Vector{T}(undef, nψ)
    return MoCSolution{T}(zero(T), zero(T), 0, φ, φ_prev, q, boundary_ψ, start_boundary_ψ)
end

function solve(prob::MoCProblem; max_iterations::Int=1000, max_residual::Real=1e-7)
    return _solve_eigenvalue_problem(prob, max_iterations, max_residual)
end

function _solve_eigenvalue_problem(prob::MoCProblem, max_iter::Int, max_ϵ::T) where {T<:Real}

    sol = MoCSolution{T}(prob.nφ, prob.nψ)

    optical_length!(prob)

    @set! sol.keff = one(T)
    set_uniform_φ!(sol, one(T))
    update_prev_φ!(sol)
    set_uniform_start_boundary_ψ!(sol, zero(T))
    update_boundary_ψ!(sol)

    @set! sol.iterations = 0
    while sol.iterations < max_iter

        normalize_fluxes!(sol, prob)
        compute_q!(sol, prob)
        compute_φ!(sol, prob)
        @set! sol.keff *= total_fission_source(sol, prob)
        @set! sol.residual = residual(sol, prob)
        update_prev_φ!(sol)

        # we need at least three iterations:
        #  0. boundary conditions do not exists (unless everything is vaccum)
        #  1. we can compute a residual but it is not correct (previous iteration is fruit)
        #  2. now we can compute a correct residual

        if sol.iterations > 1 && isless(sol.residual, max_ϵ)
            break
        end

        @show sol.residual

        @set! sol.iterations += 1
    end

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
    @unpack materials, cell_tag, tag_to_idx = prob

    for segment in track.segments
        @unpack ℓ, τ = segment

        resize!(τ, NGroups)

        #! no creo que element sea igual al id bien que va con 1:nfsr
        i = segment.element

        # dada una cell id, encuentro el tag de Gridap
        tag = cell_tag[i]

        # dado un tag de Gridap, busco el idx de NeutronTransport
        idx = tag_to_idx[tag]

        # tomo el material
        material = materials[idx]

        @unpack Σt = material

        for g in 1:NGroups
            τ[g] = Σt[g] * ℓ
        end
    end
end

@inline set_uniform_φ!(sol::MoCSolution, φ::Real) = fill!(sol.φ, φ)
@inline set_uniform_start_boundary_ψ!(sol::MoCSolution, ψ::Real) = fill!(sol.start_boundary_ψ, ψ)
@inline update_prev_φ!(sol::MoCSolution) = copy!(sol.φ_prev, sol.φ)
@inline update_boundary_ψ!(sol::MoCSolution) = copy!(sol.boundary_ψ, sol.start_boundary_ψ)

function normalize_fluxes!(sol::MoCSolution, prob::MoCProblem)
    @unpack φ, φ_prev, boundary_ψ, start_boundary_ψ = sol

    # total fission source
    qft = total_fission_source(sol, prob)

    # λ is the normalization factor
    λ = 1 / qft
    φ .*= λ
    φ_prev .*= λ
    boundary_ψ .*= λ
    start_boundary_ψ .*= λ

    return nothing
end

function total_fission_source(sol::MoCSolution{T}, prob::MoCProblem) where {T}
    NGroups = ngroups(prob)
    NRegions = nregions(prob)
    @unpack φ = sol
    @unpack trackgenerator, materials, cell_tag, tag_to_idx = prob
    @unpack volumes = trackgenerator

    qft = zero(T)
    for i in 1:NRegions

        # dada una cell id, encuentro el tag de Gridap
        tag = cell_tag[i]

        # dado un tag de Gridap, busco el idx de NeutronTransport
        idx = tag_to_idx[tag]

        # tomo el material
        material = materials[idx]

        @unpack νΣf = material

        for g′ in 1:NGroups
            ig′ = @region_index(i, g′)
            qft += νΣf[g′] * φ[ig′] * volumes[i]
        end
        # qft += sum(
        #     νΣf[g′] * φ[@region_index(i, g′)] * volumes[i] for g′ in 1:NGroups
        # )
    end

    return qft
end

function compute_q!(sol::MoCSolution{T}, prob::MoCProblem) where {T}
    NGroups = ngroups(prob)
    NRegions = nregions(prob)
    @unpack keff, φ, q = sol
    @unpack materials, cell_tag, tag_to_idx = prob

    for i in 1:NRegions

        # dada una cell id `i`, encuentro el tag de Gridap
        tag = cell_tag[i]

        # dado un tag de Gridap, busco el idx de NeutronTransport
        idx = tag_to_idx[tag]

        # tomo el material
        material = materials[idx]

        @unpack χ, Σt, νΣf, Σs0 = material

        for g in 1:NGroups
            ig = @region_index(i, g)
            qig = zero(T)
            for g′ in 1:NGroups
                ig′ = @region_index(i, g′)
                qig += Σs0[g′, g] * φ[ig′]
                qig += 1 / keff * χ[g] * νΣf[g′] * φ[ig′]
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
    @unpack φ, q = sol
    @unpack quadrature = prob

    n_polar_2 = npolar2(quadrature.polar)

    t = universal_id(track)
    d = Int32(dir)

    i = @angular_index(t, d, 1, 1)
    j = i + NGroups * n_polar_2 - 1
    boundary_ψ = @view sol.boundary_ψ[i:j]

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

    #! no creo que element sea igual al id bien que va con 1:nfsr
    i = segment.element
    a = track.azim_idx
    n_polar_2 = npolar2(polar)

    # TODO: es este el mejor orden posible?
    for p in 1:n_polar_2, g in 1:NGroups
        pg = @reduced_angular_index(p, g)
        ig = @region_index(i, g)
        Δψ = (boundary_ψ[pg] - q[ig]) * (1 - exp(-τ[g] / sinθs[p]))
        φ[ig] += 2 * ω[a, p] * Δψ
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
        flag = !isequal(boundary_out(next_track), Vaccum)
    elseif dir == Backward
        next_track = current_track.next_track_bwd
        next_track_dir = dir_next_track_bwd(current_track)
        flag = !isequal(boundary_in(next_track), Vaccum)
    end

    t = universal_id(next_track)
    d = Int32(next_track_dir)

    # TODO: es este el mejor orden posible?
    for p in 1:n_polar_2, g in 1:NGroups
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
    @unpack trackgenerator, materials, cell_tag, tag_to_idx = prob
    @unpack volumes = trackgenerator

    for i in 1:NRegions

        # dada una cell id `i`, encuentro el tag de Gridap
        tag = cell_tag[i]

        # dado un tag de Gridap, busco el idx de NeutronTransport
        idx = tag_to_idx[tag]

        # tomo el material
        material = materials[idx]

        @unpack Σt = material

        for g in 1:NGroups
            ig = @region_index(i, g)
            φ[ig] /= (Σt[g] * volumes[i])
            φ[ig] += 4π * q[ig]
        end
    end

    return nothing
end

function residual(sol::MoCSolution{T}, prob::MoCProblem) where {T}
    NGroups = ngroups(prob)
    NRegions = nregions(prob)
    @unpack keff, φ, φ_prev = sol
    @unpack materials, cell_tag, tag_to_idx = prob

    # calculamos la fuente total (integrada en energia g) para cada celda i (fuente new y
    # old). Podriamos almacenar la new y dejarla como old para la siguiente pasada...

    ϵ = zero(T)
    for i in 1:NRegions

        # dada una cell id, encuentro el tag de Gridap
        tag = cell_tag[i]

        # dado un tag de Gridap, busco el idx de NeutronTransport
        idx = tag_to_idx[tag]

        # tomo el material
        material = materials[idx]

        new_qi = zero(T)
        old_qi = zero(T)

        @unpack νΣf, Σs0 = material

        # total fission source in each region (χ sum 1 when ∫ in g)
        for g′ in 1:NGroups
            ig′ = @region_index(i, g′)
            νΣfg′ = νΣf[g′]
            new_qi += νΣfg′ * φ[ig′]
            old_qi += νΣfg′ * φ_prev[ig′]
        end
        # new_qi = sum(νΣf[g′] * φ[ig′] for g′ in 1:NGroups)
        # old_qi = sum(νΣf[g′] * φ_pre[ig′] for g′ in 1:NGroups)

        new_qi /= keff
        old_qi /= keff

        # total scattering source
        for g in 1:NGroups, g′ in 1:NGroups
            ig′ = @region_index(i, g′)
            Σsgg′ = Σs0[g′, g]
            old_qi += Σsgg′ * φ[ig′]
            new_qi += Σsgg′ * φ_prev[ig′]
        end

        if old_qi > 0
            ϵ += ((new_qi - old_qi) / old_qi)^2
        end
    end

    ϵ = sqrt(ϵ / NRegions)

    return ϵ
end
