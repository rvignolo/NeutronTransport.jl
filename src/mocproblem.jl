
struct MoCProblem{T<:Real,G<:TrackGenerator,Q<:Quadrature,M}
    # as part of the type?
    dimension::Int
    groups::Int
    nfsr::Int

    max_iterations::Int
    max_residual::T

    # scalar flux
    nφ::Int
    φ::Vector{T}
    φ_prev::Vector{T}

    q::Vector{T}  # reduced source

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
    tg::TrackGenerator{M,Q,T1}, pq::PolarQuadrature{N,T2}, materials::P;
    max_iter::Int=1000, max_residual::Real=1e-7
) where {M,Q,T1,N,T2,P}
    @unpack mesh, azimuthal_quadrature, n_total_tracks = tg

    dimension = num_dims(mesh)
    groups = 2 # from materials? sip.
    nfsr = num_cells(mesh)
    nφ = nfsr * groups
    nψ = n_total_tracks * npolar(pq) * groups

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
        gridap_tag = get_tag_from_name(key)
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
        (i - 1) * groups + g
    end
    return esc(ex)
end
# puedo hacerlas como function con @inline pero necesitan `groups` o `nfsr`

macro angular_index(t, d, p, g)
    ex = quote
        (t - 1) * 2 * n_polar_2 * groups + (d - 1) * n_polar_2 * groups + (p - 1) * groups + g
    end
    return esc(ex)
end

macro reduced_angular_index(p, g)
    ex = quote
        (p - 1) * groups + g
    end
    return esc(ex)
end

solve(moc::MoCProblem) = _solve_eigenvalue_problem(moc)

function _solve_eigenvalue_problem(moc::MoCProblem{T}) where {T<:Real}

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

        if iter > 1 && isless(residual, moc.max_residual)
            break
        end

        iter += 1
    end

end

@inline set_uniform_φ!(moc::MoCProblem, φ::Real) = fill!(moc.φ, φ)
@inline set_uniform_start_boundary_ψ!(moc::MoCProblem, ψ::Real) = fill!(moc.start_boundary_ψ, ψ)
@inline update_prev_φ!(moc::MoCProblem) = copy!(moc.φ_prev, moc.φ)
@inline update_boundary_ψ!(moc::MoCProblem) = copy!(moc.boundary_ψ, moc.start_boundary_ψ)

function normalize_fluxes!(moc::MoCProblem{T}) where {T}
    @unpack groups, nfsr, φ, materials, cell_tag, tag_to_idx = moc

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

function total_fission_source(moc::MoCProblem)
    @unpack groups, nfsr, φ, materials, cell_tag, tag_to_idx = moc

    qft = zero(T)
    for i in 1:nfsr

        # dada una cell id, encuentro el tag de Gridap
        tag = cell_tag[i]

        # dado un tag de Gridap, busco el idx de NeutronTransport
        idx = tag_to_idx[tag]

        # tomo el material
        material = materials[idx]

        #! lo tengo que crear, quizas calculando los volumes con los tracks... o no
        volume = cell_volume[i]

        @unpack νΣf = material

        for g′ in 1:groups
            ig′ = @region_index(i, g′)
            qft += νΣf[g′] * φ[ig′] * volume
        end
        # qft += sum(
        #     νΣf[g′] * φ[@region_index(i, g′)] * volume for g′ in 1:groups
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

            qig /= (4 * π * Σt[g])
            q[ig] = qig
        end
    end

    return nothing
end



function compute_φ!(moc::MoCProblem{T}) where {T}

    set_uniform_φ!(moc, zero(T))
    update_boundary_ψ!(moc)

    for (t, track) in enumerate(trackgenerator.tracks)
        a = @angular_index(t, 0, 0, 0)
        b = a + groups * n_polar_2
        boundary_ψ = @view moc.boundary_ψ[a:b]
        for segment in track.segments
            tally_φ!(moc, segment, boundary_ψ, i) # en lugar de i mandar track y ese tiene el uid que es t y el azim index, entonces el @view lo hago en tally? ah pero abajo lo necesito de nuevo al @view en set_start_boundary_ψ
        end

        # TODO: transferir los flujos angulares fwd al track correspondiente
        set_start_boundary_ψ(moc, )

        a = @angular_index(t, 1, 0, 0)
        b = a + groups * n_polar_2
        boundary_ψ = @view moc.boundary_ψ[a:b]
        for segment in reverse!(track.segments)
            tally_φ!(moc, segment, boundary_ψ, i) # en lugar de i mandar track y ese tiene el uid que es t y el azim index, entonces el @view lo hago en tally? ah pero abajo lo necesito de nuevo al @view
        end

        reverse!(track.segments)

    end

end


function tally_φ!(moc, segment, boundary_ψ, azim_idx)
    for g in 1:groups
        for p in 1:n_polar_2
            exponential = 1. # TODO
            pg = @reduced_angular_index(p, g)
            ig = @region_index(segment.element, g) # element deberia llamarse cell_id?
            Δψ = (boundary_ψ[pg] - moc.q[ig]) * exponential
            φ[ig] += 2 * moc.quadrature.ω[azim_idx, p] * Δψ
            boundary_ψ[pg] -= Δψ
        end
    end
    return nothing
end

multiplication_factor(moc::MoCProblem, keff::Real) = total_fission_source(moc) * keff

function residual(moc::MoCProblem{T}, keff::Real) where {T}
    @unpack nfsr, groups, φ, φ_prev = moc

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