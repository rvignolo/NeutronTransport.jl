
struct MoCProblem{T<:Real,G<:TrackGenerator,Q<:Quadrature,M}
    # as part of the type?
    dimension::Int
    groups::Int
    nfsr::Int

    max_iterations::Int
    max_residual::T

    # scalar flux
    nφ::Int
    φ::Vector{T} # me parece que algo de Gridap va a ser mejor
    φ_prev::Vector{T}

    reduced_souce::Vector{T}

    # angular flux
    nψ::Int
    boundary_ψ::Vector{T}
    start_boundary_ψ::Vector{T}

    trackgenerator::G
    quadrature::Q

    materials::M
    cells_tag::Vector{Int8}  # gridap cell to gridap tag
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
    reduced_source = Vector{T}(undef, nφ)
    boundary_ψ = Vector{T}(undef, nψ)
    start_boundary_ψ = Vector{T}(undef, nψ)

    #! TODO: mapear dado un indice de celda y un indice de grupo, el φ index que le corresponde

    quadrature = Quadrature(azimuthal_quadrature, pq)

    face_labeling = get_face_labeling(tg.mesh.model)
    cells_tag = get_face_tag(face_labeling, dimension)
    tag_to_name = face_labeling.tag_to_name
    tag_to_idx = Vector{Int8}(undef, length(tag_to_name))

    # me van a quedar sin definir las posiciones q no correspondan a materiales...
    for (i, pair) in enumerate(materials)
        key, value = pair
        # dado un name mio, busco su idx en gridap
        gridap_idx = findfirst(x -> isequal(x, key), tag_to_name)
        tag_to_idx[gridap_idx] = i
    end

    # mmm. no necesito que sea named tuple ya que todos los types son iguales y puedo ir a
    # un vector. por otro lado, tambien podria usar una tupla, no es necesario la named part
    # materials_nt = (; (Symbol(k) => v for (k, v) in materials)...)
    materials_v = [values(materials)...]

    return MoCProblem(
        dimension, groups, nfsr, max_iter, max_residual, nφ, φ, φ_prev, reduced_source, nψ,
        boundary_ψ, start_boundary_ψ, tg, quadrature, materials_v, cells_tag, tag_to_idx
    )
end

macro angular_index(t, d, p, g)
    return :(t * 2 * n_polar_2 * groups + d * n_polar_2 * groups + p * groups + g)
end

solve(moc::MoCProblem) = _solve_eigenvalue_problem(moc)

function _solve_eigenvalue_problem(moc::MoCProblem{T}) where {T<:Real}

    #! donde lo almacenamos? seguro tienen que ser visible para todos
    keff = one(T)

    set_uniform_φ!(moc, one(T))
    update_prev_φ!(moc)
    set_uniform_start_boundary_ψ!(moc, zero(T))
    update_boundary_ψ!(moc)

    iter = 0
    while iter < moc.max_iterations

        normalize_fluxes!(moc)
        compute_q!(moc)
        compute_φ!(moc)
        keff = multiplication_factor(moc) #! ojo este se debe usar en calculos!
        ϵ = residual(moc)
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

    @unpack mesh = moc.trackgenerator

    total_fission_source = zero(T)

    for i in 1:moc.nfsr

        #! no existe esto aun
        #! dada una cell tenemos que conseguir el material, i.e. las secciones eficaces
        cell = mesh.cells[i]
        material = cell.material

        # TODO: check if for loop is faster
        total_fission_source += sum(
            material.xs.νΣf[g′] * φ[cell.index[g′]] * cell.volume for g′ in 1:moc.groups
        )

    end

    normalization_factor = 1 / total_fission_source
    φ .*= normalization_factor
    φ_prev .*= normalization_factor
    boundary_ψ .*= normalization_factor
    start_boundary_ψ .*= normalization_factor

    return nothing
end

function compute_q!(moc::MoCProblem{T}) where {T}
    @unpack groups, materials, cells_tag, tag_to_idx = moc

    for i in 1:moc.nfsr

        # dada una cell id, encuentro el tag de Gridap
        tag = moc.cells_tag[i]

        # dado un tag de Gridap, busco el idx de NeutronTransport
        idx = moc.tag_to_idx[tag]

        # tomo el material
        material = moc.materials[idx]

        @unpack χ, Σt, νΣf, Σs0 = material

        for g in 1:groups

            #! use view instead
            reduced_source[cell.index[g]] = zero(T)

            for g′ in 1:groups
                reduced_source[cell.index[g]] += Σs0[g′][g] * φ[cell.index[g_prime]]
                reduced_source[cell.index[g]] += 1 / keff * χ[g] * material.xs.νΣf[g′] * φ[cell.index[g_prime]]
            end

            #! aca se computa Σtotal, en realidad ya lo deberia hacer en el constructor

            reduced_source[cell.index[g]] /= (4 * π * Σt[g])
        end

    end

    return nothing
end

function compute_φ(moc::MoCProblem) end
