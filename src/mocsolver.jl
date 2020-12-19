
struct MoCSolver{T<:Real,G<:TrackGenerator,PQ<:PolarQuadrature,Q<:Quadrature}
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
    polar_quadrature::PQ
    quadrature::Q
end

function MoCSolver(
    tg::TrackGenerator{M,Q,T1}, pq::PolarQuadrature{N,T2};
    max_iter::Int=1000, max_residual::Real=1e-7
) where {M,Q,T1,N,T2}
    @unpack mesh, n_total_tracks = tg

    dimension = num_dims(mesh)
    groups = 2 # from materials?
    nfsr = num_cells(mesh)
    nφ = nfsr * groups
    nψ = n_total_tracks * npolar(pq) * groups

    T = promote_type(T1, T2)
    φ = Vector{T}(undef, nφ)
    φ_prev = Vector{T}(undef, nφ)
    reduced_source = Vector{T}(undef, nφ)
    boundary_ψ = Vector{T}(undef, nψ)
    start_boundary_ψ = Vector{T}(undef, nψ)

    quadrature = Quadrature(tg, pq)

    return MoCSolver(
        dimension, groups, nfsr, max_iter, max_residual, nφ, φ, φ_prev,
        reduced_source, nψ, boundary_ψ, start_boundary_ψ, tg, pq, quadrature
    )
end

macro angular_index(t, d, p, g)
    return :(t * 2 * n_polar_2 * groups + d * n_polar_2 * groups + p * groups + g)
end

solve(moc::MoCSolver) = _solve_eigenvalue_problem(moc)

function _solve_eigenvalue_problem(moc::MoCSolver{T}) where {T<:Real}

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

@inline set_uniform_φ!(moc::MoCSolver, φ::Real) = fill!(moc.φ, φ)
@inline set_uniform_start_boundary_ψ!(moc::MoCSolver, ψ::Real) = fill!(moc.start_boundary_ψ, ψ)
@inline update_prev_φ!(moc::MoCSolver) = copy!(moc.φ_prev, moc.φ)
@inline update_boundary_ψ!(moc::MoCSolver) = copy!(moc.boundary_ψ, moc.start_boundary_ψ)

function normalize_fluxes!(moc::MoCSolver{T}) where {T}

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

function compute_q!(moc::MoCSolver{T}) where {T}

    for i in 1:moc.nfsr

        #! no existe esto aun
        #! dada una cell tenemos que conseguir el material, i.e. las secciones eficaces
        cell = mesh.cells[i]
        material = cell.material

        for g in 1:groups

            # use view instead

            reduced_source[cell.index[g]] = zero(T)

            for g′ in 1:groups
                reduced_source[cell.index[g]] += material.xs.νΣS0[g′][g] * φ[cell.index[g_prime]]
                reduced_source[cell.index[g]] += 1 / keff * χ[g] * material.xs.νΣf[g′] * φ[cell.index[g_prime]]
            end

            #! aca se computa Σtotal


            reduced_source[cell.index[g]] /= (4 * π * material.xs.Σt[g])
        end

    end

    return nothing
end

function compute_φ(moc::MoCSolver) end
