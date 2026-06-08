@testset "MoC smoke solve with local optical lengths" begin
    cases = (
        ("vacuum", BoundaryConditions()),
        ("reflective", BoundaryConditions(
            top=Reflective, bottom=Reflective, left=Reflective, right=Reflective
        )),
    )

    for (name, bcs) in cases
        @testset "$name boundaries" begin
            prob, sol, tg = solve_single_cell(bcs)

            @test isfinite(sol.keff)
            @test isfinite(sol.residual)
            @test sol.iterations == 3
            @test length(prob.optical_lengths) == tg.n_total_tracks
            @test length(prob.attenuation_factors) == tg.n_total_tracks
            @test all(τ -> size(τ, 1) == 2, prob.optical_lengths)
            @test all(τ -> all(isfinite, τ), prob.optical_lengths)
            @test all(A -> size(A, 1) == 2, prob.attenuation_factors)
            @test all(A -> size(A, 2) == 1, prob.attenuation_factors)
            @test all(A -> all(isfinite, A), prob.attenuation_factors)
            @test all(A -> all(a -> zero(a) <= a <= one(a), A), prob.attenuation_factors)

            for (uid, track) in enumerate(tg.tracks_by_uid)
                @test size(prob.optical_lengths[uid], 2) == length(track.segments)
                @test size(prob.attenuation_factors[uid], 3) == length(track.segments)
                expected = @. -expm1(
                    -prob.optical_lengths[uid] / prob.quadrature.polar.sinθs[1]
                )
                @test prob.attenuation_factors[uid][:, 1, :] ≈ expected
            end
        end
    end
end

@testset "core MoC kernels infer concrete return types" begin
    import RayTracing: Forward, universal_id

    prob, tg = single_cell_problem(BoundaryConditions(
        top=Reflective, bottom=Reflective, left=Reflective, right=Reflective
    ))
    sol = NeutronTransport.MoCSolution{eltype(prob)}(prob)
    NeutronTransport.optical_length!(prob)
    NeutronTransport.set_uniform_φ!(sol, one(eltype(prob)))
    NeutronTransport.set_uniform_Q!(sol, zero(eltype(prob)))
    NeutronTransport.set_uniform_start_boundary_ψ!(sol, zero(eltype(prob)))
    NeutronTransport.update_boundary_ψ!(sol)
    NeutronTransport.normalize_fluxes!(sol, prob)
    NeutronTransport.compute_q!(sol, prob)

    track = tg.tracks_by_uid[1]
    segment = track.segments[1]
    t = universal_id(track)
    attenuation = prob.attenuation_factors[t]
    n_polar_half_count = NeutronTransport.n_polar_half(prob.quadrature.polar)
    NGroups = NeutronTransport.ngroups(prob)
    boundary_offset = (t - 1) * 2 * n_polar_half_count * NGroups +
        Int(Forward) * n_polar_half_count * NGroups

    @test (@inferred NeutronTransport.compute_q_region!(
        sol.q, sol.φ, prob, 1, inv(sol.keff)
    )) === nothing
    @test (@inferred NeutronTransport.tally_φ!(
        sol, prob, track, segment, attenuation, 1,
        sol.boundary_ψ, boundary_offset, sol.φ
    )) === nothing
    @test (@inferred NeutronTransport.tally!(sol, prob, track, Forward)) === nothing
    @test (@inferred NeutronTransport.set_start_boundary_ψ!(
        sol, prob, track, sol.boundary_ψ, boundary_offset, Forward
    )) === nothing
    @test (@inferred NeutronTransport.total_fission_source(sol, prob)) isa eltype(prob)
    @test (@inferred NeutronTransport.compute_q!(sol, prob)) === nothing
    @test (@inferred NeutronTransport.compute_φ!(sol, prob)) === nothing
    @test (@inferred NeutronTransport.residual(sol, prob)) isa eltype(prob)

    if Base.Threads.nthreads() > 1
        @test (@inferred NeutronTransport.total_fission_source(
            sol, prob; parallel=true
        )) isa eltype(prob)
        @test (@inferred NeutronTransport.compute_q!(sol, prob; parallel=true)) === nothing
        @test (@inferred NeutronTransport.compute_φ!(sol, prob; parallel=true)) === nothing
        @test (@inferred NeutronTransport.residual(sol, prob; parallel=true)) isa eltype(prob)
    end
end

@testset "region-parallel solve matches serial solve" begin
    model = openmoc_reflective_grid_model()
    materials = [openmoc_homogeneous_medium(), openmoc_water_medium()]

    serial_prob = reflected_problem(model, materials)
    parallel_prob = reflected_problem(model, materials)

    @test NeutronTransport.validate_parallel_sweep_boundaries(parallel_prob) === nothing

    serial_sol = solve(serial_prob; max_iterations=5, max_residual=0.0)
    parallel_sol = solve(parallel_prob; max_iterations=5, max_residual=0.0, parallel=true)

    @test parallel_sol.keff ≈ serial_sol.keff
    @test parallel_sol.residual ≈ serial_sol.residual
    @test parallel_sol.φ ≈ serial_sol.φ
    @test parallel_sol.Q ≈ serial_sol.Q
end

@testset "material-coalesced FSR mapping" begin
    materials = [openmoc_homogeneous_medium(), openmoc_water_medium()]
    bcs = BoundaryConditions(
        top=Reflective, bottom=Reflective, left=Reflective, right=Reflective
    )
    tg = TrackGenerator(openmoc_reflective_grid_model(), 8, 0.25;
        bcs=bcs, volume_correction=true
    )
    trace!(tg)
    segmentize!(tg)

    default_prob = MoCProblem(tg, TabuchiYamamoto(2), materials)
    cell_to_fsr = cell_material_ids(tg, materials)
    coalesced_prob = MoCProblem(
        tg, TabuchiYamamoto(2), materials; cell_to_fsr=cell_to_fsr
    )

    @test default_prob.cell_to_fsr == Int32.(1:9)
    @test default_prob.fsr_tag == cell_to_fsr
    @test default_prob.volumes ≈ tg.volumes

    @test length(coalesced_prob.fsr_tag) == 2
    @test coalesced_prob.fsr_tag == Int32[1, 2]
    @test count(==(1), coalesced_prob.cell_to_fsr) == 1
    @test count(==(2), coalesced_prob.cell_to_fsr) == 8
    @test coalesced_prob.volumes[1] ≈ sum(tg.volumes[cell_to_fsr .== 1])
    @test coalesced_prob.volumes[2] ≈ sum(tg.volumes[cell_to_fsr .== 2])

    sol = solve(coalesced_prob; max_iterations=2, max_residual=0.0)
    @test length(sol.φ) == 2 * 2
    @test isfinite(sol.keff)
    @test all(isfinite, sol.φ)
end

@testset "invalid FSR mappings" begin
    materials = [openmoc_homogeneous_medium(), openmoc_water_medium()]
    bcs = BoundaryConditions(
        top=Reflective, bottom=Reflective, left=Reflective, right=Reflective
    )
    tg = TrackGenerator(openmoc_reflective_grid_model(), 8, 0.25;
        bcs=bcs, volume_correction=true
    )
    trace!(tg)
    segmentize!(tg)

    @test_throws ArgumentError MoCProblem(
        tg, TabuchiYamamoto(2), materials; cell_to_fsr=fill(1, 8)
    )
    @test_throws ArgumentError MoCProblem(
        tg, TabuchiYamamoto(2), materials; cell_to_fsr=[1, fill(3, 8)...]
    )
    @test_throws ArgumentError MoCProblem(
        tg, TabuchiYamamoto(2), materials; cell_to_fsr=fill(1, 9)
    )

    prob = MoCProblem(
        tg, TabuchiYamamoto(2), materials; cell_to_fsr=fill(1, 9), fsr_to_xs=[1]
    )
    @test length(prob.fsr_tag) == 1
    @test prob.volumes[1] ≈ sum(tg.volumes)
end

@testset "zero-volume FSR validation" begin
    prob, tg = single_cell_problem(BoundaryConditions())
    prob.volumes[1] = 0.0

    @test_throws ArgumentError NeutronTransport.solve(prob; max_iterations=1)
end
