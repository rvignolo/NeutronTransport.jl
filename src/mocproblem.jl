
struct MoCProblem{T<:Real,G<:TrackGenerator,Q<:Quadrature,M}
    # TODO: as part of the type? si
    dimension::Int
    groups::Int
    nfsr::Int

    max_iterations::Int
    max_residual::T

    # scalar flux
    nφ::Int
    φ::Vector{T}
    φ_prev::Vector{T}

    # reduced source
    q::Vector{T}

    # angular flux
    nψ::Int
    boundary_ψ::Vector{T}
    start_boundary_ψ::Vector{T}

    trackgenerator::G
    quadrature::Q

    materials::M
    cell_tag::Vector{Int8}  # gridap cell to gridap tag
    tag_to_idx::Vector{Int8} #  gridap tag to NeutronTransport material index
end

# puedo recibir dict y pasarlo aca a tupla usando las mismas tags que gridap (incluso como
# los materiales tienen las mismas dimensiones, puedo usar vectores de structs y tener un
# type defindo...). Ver si NamedTuples.jl tiene el constructor de named tuple a partir de dict.
function MoCProblem(
    tg::TrackGenerator{T1}, pq::PolarQuadrature{N,T2}, materials::P;
    max_iter::Int=1000, max_residual::Real=1e-7
) where {T1,N,T2,P}
    @unpack mesh, azimuthal_quadrature, n_total_tracks = tg

    dimension = num_dims(mesh)
    groups = ngroups.(values(materials))[1] # TODO: check all equal
    nfsr = num_cells(mesh)
    nφ = nfsr * groups
    nψ = n_total_tracks * 2 * npolar2(pq) * groups

    T = promote_type(T1, T2)
    φ = Vector{T}(undef, nφ)
    φ_prev = Vector{T}(undef, nφ)
    q = Vector{T}(undef, nφ)
    boundary_ψ = Vector{T}(undef, nψ)
    start_boundary_ψ = Vector{T}(undef, nψ)

    quadrature = Quadrature(azimuthal_quadrature, pq)

    face_labeling = get_face_labeling(tg.mesh.model)
    cell_tag = get_face_tag(face_labeling, dimension)
    tag_to_name = face_labeling.tag_to_name
    tag_to_idx = Vector{Int8}(undef, length(tag_to_name))

    # me van a quedar sin definir las posiciones q no correspondan a materiales...
    for (i, pair) in enumerate(materials)
        key, value = pair
        # dado un name mio, busco su tag en gridap
        gridap_tag = get_tag_from_name(face_labeling, key)
        tag_to_idx[gridap_tag] = i
    end

    # mmm. no necesito que sea named tuple ya que todos los types son iguales y puedo ir a
    # un vector. por otro lado, tambien podria usar una tupla, no es necesario la named part
    # materials_nt = (; (Symbol(k) => v for (k, v) in materials)...)
    materials_v = [values(materials)...]

    return MoCProblem(
        dimension, groups, nfsr, max_iter, max_residual, nφ, φ, φ_prev, q, nψ,
        boundary_ψ, start_boundary_ψ, tg, quadrature, materials_v, cell_tag, tag_to_idx
    )
end

# si quiero que los primeros nfsr sean para el 1er group, etc...
# macro region_index(i, g)
#     return :((g - 1) * nfsr + i)
# end

# si quiero que esten intercalados
macro region_index(i, g)
    ex = quote
        ($i - 1) * groups + $g
    end
    return esc(ex)
end
# puedo hacerlas como function con @inline pero necesitan `groups` o `nfsr`

macro angular_index(t, d, p, g)
    ex = quote
        ($t - 1) * 2 * n_polar_2 * groups + $d * n_polar_2 * groups + ($p - 1) * groups + $g
    end
    return esc(ex)
end

macro reduced_angular_index(p, g)
    ex = quote
        ($p - 1) * groups + $g
    end
    return esc(ex)
end

solve(moc::MoCProblem) = _solve_eigenvalue_problem(moc)

function _solve_eigenvalue_problem(moc::MoCProblem{T}) where {T<:Real}

    optical_length!(moc)

    keff = one(T)
    set_uniform_φ!(moc, one(T))
    update_prev_φ!(moc)
    set_uniform_start_boundary_ψ!(moc, zero(T))
    update_boundary_ψ!(moc)

    iter = 0
    while iter < moc.max_iterations

        normalize_fluxes!(moc)
        compute_q!(moc, keff)
        compute_φ!(moc)
        keff = multiplication_factor(moc, keff)
        ϵ = residual(moc, keff)
        update_prev_φ!(moc)

        # we need at least three iterations:
        #  0. boundary conditions do not exists (unless everything is vaccum)
        #  1. we can compute a residual but it is not correct (previous iteration is fruit)
        #  2. now we can compute a correct residual

        if iter > 1 && isless(ϵ, moc.max_residual)
            break
        end

        @show ϵ

        iter += 1
    end

end

function optical_length!(moc)
    @unpack trackgenerator = moc
    @unpack tracks_by_uid = trackgenerator

    for track in tracks_by_uid
        _optical_length!(moc, track)
    end

    return nothing
end

function _optical_length!(moc, track)
    @unpack groups, materials, cell_tag, tag_to_idx = moc

    for segment in track.segments
        @unpack ℓ, τ = segment
        resize!(τ, groups)

        #! no creo que element sea igual al id bien que va con 1:nfsr
        i = segment.element

        # dada una cell id, encuentro el tag de Gridap
        tag = cell_tag[i]

        # dado un tag de Gridap, busco el idx de NeutronTransport
        idx = tag_to_idx[tag]

        # tomo el material
        material = materials[idx]

        @unpack Σt = material

        for g in 1:groups
            τ[g] = Σt[g] * ℓ
        end
    end
end

@inline set_uniform_φ!(moc::MoCProblem, φ::Real) = fill!(moc.φ, φ)
@inline set_uniform_start_boundary_ψ!(moc::MoCProblem, ψ::Real) = fill!(moc.start_boundary_ψ, ψ)
@inline update_prev_φ!(moc::MoCProblem) = copy!(moc.φ_prev, moc.φ)
@inline update_boundary_ψ!(moc::MoCProblem) = copy!(moc.boundary_ψ, moc.start_boundary_ψ)

function normalize_fluxes!(moc::MoCProblem)
    @unpack φ, φ_prev, boundary_ψ, start_boundary_ψ = moc

    # total fission source
    qft = total_fission_source(moc)

    # λ is the normalization factor
    λ = 1 / qft
    φ .*= λ
    φ_prev .*= λ
    boundary_ψ .*= λ
    start_boundary_ψ .*= λ

    return nothing
end

function total_fission_source(moc::MoCProblem{T}) where {T}
    @unpack groups, nfsr, φ, trackgenerator, materials, cell_tag, tag_to_idx = moc
    @unpack volumes = trackgenerator

    qft = zero(T)
    for i in 1:nfsr

        # dada una cell id, encuentro el tag de Gridap
        tag = cell_tag[i]

        # dado un tag de Gridap, busco el idx de NeutronTransport
        idx = tag_to_idx[tag]

        # tomo el material
        material = materials[idx]

        @unpack νΣf = material

        for g′ in 1:groups
            ig′ = @region_index(i, g′)
            qft += νΣf[g′] * φ[ig′] * volumes[i]
        end
        # qft += sum(
        #     νΣf[g′] * φ[@region_index(i, g′)] * volumes[i] for g′ in 1:groups
        # )
    end

    return qft
end

function compute_q!(moc::MoCProblem{T}, keff::Real) where {T}
    @unpack groups, nfsr, φ, q, materials, cell_tag, tag_to_idx = moc

    for i in 1:nfsr

        # dada una cell id `i`, encuentro el tag de Gridap
        tag = cell_tag[i]

        # dado un tag de Gridap, busco el idx de NeutronTransport
        idx = tag_to_idx[tag]

        # tomo el material
        material = materials[idx]

        @unpack χ, Σt, νΣf, Σs0 = material

        for g in 1:groups
            ig = @region_index(i, g)
            qig = zero(T)
            for g′ in 1:groups
                ig′ = @region_index(i, g′)
                qig += Σs0[g′, g] * φ[ig′]
                qig += 1 / keff * χ[g] * νΣf[g′] * φ[ig′]
            end
            # qig = sum(
            #     Σs0[g′, g] * φ[@region_index(i, g′)] +
            #     1 / keff * χ[g] * νΣf[g′] * φ[@region_index(i, g′)]
            #     for g′ in 1:groups
            # )

            qig /= (4π * Σt[g])
            q[ig] = qig
        end
    end

    return nothing
end

function compute_φ!(moc::MoCProblem{T}) where {T}
    @unpack groups, φ, q, trackgenerator, quadrature = moc
    @unpack tracks_by_uid = trackgenerator

    n_polar_2 = npolar2(quadrature.polar)

    set_uniform_φ!(moc, zero(T))
    update_boundary_ψ!(moc)

    for track in tracks_by_uid
        tally(moc, track, RayTracing.Forward)
        tally(moc, track, RayTracing.Backward)
    end

    add_q_to_φ!(moc)

    return nothing
end

function tally(moc, track, dir)
    @unpack groups, φ, q, quadrature = moc

    n_polar_2 = npolar2(quadrature.polar)

    t = RayTracing.universal_id(track)
    d = Int32(dir)

    i = @angular_index(t, d, 1, 1)
    j = i + groups * n_polar_2 - 1
    boundary_ψ = @view moc.boundary_ψ[i:j]

    segments = dir == RayTracing.Forward ? track.segments : reverse!(track.segments)

    for segment in segments
        tally_φ!(moc, track, segment, boundary_ψ)
    end

    set_start_boundary_ψ!(moc, track, boundary_ψ, dir)

    if dir == RayTracing.Backward
        reverse!(track.segments)
    end

    return nothing
end

function tally_φ!(moc, track, segment, boundary_ψ)
    @unpack groups, φ, q, quadrature = moc
    @unpack polar, ω = quadrature
    @unpack sinθs = polar
    @unpack τ = segment

    #! no creo que element sea igual al id bien que va con 1:nfsr
    i = segment.element
    a = track.azim_idx
    n_polar_2 = npolar2(polar)

    # TODO: es este el mejor orden posible?
    for p in 1:n_polar_2, g in 1:groups
        pg = @reduced_angular_index(p, g)
        ig = @region_index(i, g)
        Δψ = (boundary_ψ[pg] - q[ig]) * (1 - exp(-τ[g] / sinθs[p]))
        φ[ig] += 2 * ω[a, p] * Δψ
        boundary_ψ[pg] -= Δψ
    end

    return nothing
end

function set_start_boundary_ψ!(moc::MoCProblem{T}, current_track, boundary_ψ, dir) where {T}
    @unpack groups, start_boundary_ψ, quadrature = moc

    n_polar_2 = npolar2(quadrature.polar)

    if dir == RayTracing.Forward
        next_track = current_track.next_track_fwd
        next_track_dir = RayTracing.dir_next_track_fwd(current_track)
        flag = !isequal(RayTracing.boundary_out(next_track), RayTracing.Vaccum)
    elseif dir == RayTracing.Backward
        next_track = current_track.next_track_bwd
        next_track_dir = RayTracing.dir_next_track_bwd(current_track)
        flag = !isequal(RayTracing.boundary_in(next_track), RayTracing.Vaccum)
    end

    t = RayTracing.universal_id(next_track)
    d = Int32(next_track_dir)

    # TODO: es este el mejor orden posible?
    for p in 1:n_polar_2, g in 1:groups
        tdpg = @angular_index(t, d, p, g)
        pg = @reduced_angular_index(p, g)
        start_boundary_ψ[tdpg] = flag ? boundary_ψ[pg] : zero(T)
    end

    return nothing
end

function add_q_to_φ!(moc)
    @unpack groups, nfsr, φ, q, trackgenerator, materials, cell_tag, tag_to_idx = moc
    @unpack volumes = trackgenerator

    for i in 1:nfsr

        # dada una cell id `i`, encuentro el tag de Gridap
        tag = cell_tag[i]

        # dado un tag de Gridap, busco el idx de NeutronTransport
        idx = tag_to_idx[tag]

        # tomo el material
        material = materials[idx]

        @unpack Σt = material

        for g in 1:groups
            ig = @region_index(i, g)
            φ[ig] /= (Σt[g] * volumes[i])
            φ[ig] += 4π * q[ig]
        end
    end

    return nothing
end

multiplication_factor(moc::MoCProblem, keff::Real) = total_fission_source(moc) * keff

function residual(moc::MoCProblem{T}, keff::Real) where {T}
    @unpack nfsr, groups, φ, φ_prev, materials, cell_tag, tag_to_idx = moc

    # calculamos la fuente total (integrada en energia g) para cada celda i (fuente new y
    # old). Podriamos almacenar la new y dejarla como old para la siguiente pasada...

    ϵ = zero(T)
    for i in 1:nfsr

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
        for g′ in 1:groups
            ig′ = @region_index(i, g′)
            νΣfg′ = νΣf[g′]
            new_qi += νΣfg′ * φ[ig′]
            old_qi += νΣfg′ * φ_prev[ig′]
        end
        # new_qi = sum(νΣf[g′] * φ[ig′] for g′ in 1:groups)
        # old_qi = sum(νΣf[g′] * φ_pre[ig′] for g′ in 1:groups)

        new_qi /= keff
        old_qi /= keff

        # total scattering source
        for g in 1:groups, g′ in 1:groups
            ig′ = @region_index(i, g′)
            Σsgg′ = Σs0[g′, g]
            old_qi += Σsgg′ * φ[ig′]
            new_qi += Σsgg′ * φ_prev[ig′]
        end

        if old_qi > 0
            ϵ += ((new_qi - old_qi) / old_qi)^2
        end
    end

    ϵ = sqrt(ϵ / nfsr)

    return ϵ
end