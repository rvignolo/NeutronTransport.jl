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
