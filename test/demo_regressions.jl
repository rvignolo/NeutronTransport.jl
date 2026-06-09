@testset "demo pincell short solve" begin
    bcs = BoundaryConditions(
        top=Reflective, bottom=Reflective, left=Reflective, right=Reflective
    )
    prob, _ = demo_problem("pincell", pincell_materials();
        n_azim=8, spacing=0.08, n_polar=2, bcs=bcs
    )
    sol = solve(prob; max_iterations=3, max_residual=0.0)

    @test count(==(1), prob.fsr_tag) > 0
    @test count(==(2), prob.fsr_tag) > 0
    @test count(==(3), prob.fsr_tag) > 0
    @test isfinite(sol.keff)
    @test sol.keff > 0
    @test isfinite(sol.residual)
    @test all(isfinite, sol.φ)
    @test sum(abs, sol.φ) > 0
end

@testset "demo BWR short solve" begin
    bcs = BoundaryConditions(
        top=Reflective, bottom=Reflective, left=Reflective, right=Reflective
    )
    prob, _ = demo_problem("bwr", bwr_materials();
        n_azim=4, spacing=0.18, n_polar=2, bcs=bcs
    )
    sol = solve(prob; max_iterations=3, max_residual=0.0)

    @test count(==(1), prob.fsr_tag) > 0
    @test count(==(2), prob.fsr_tag) > 0
    @test count(==(3), prob.fsr_tag) > 0
    @test count(==(4), prob.fsr_tag) > 0
    @test isfinite(sol.keff)
    @test sol.keff > 0
    @test isfinite(sol.residual)
    @test all(isfinite, sol.φ)
    @test sum(abs, sol.φ) > 0
end

@testset "demo C5G7 geometry short solve" begin
    bcs = BoundaryConditions(top=Vacuum, bottom=Reflective, left=Reflective, right=Vacuum)
    prob, tg = demo_problem("c5g7", c5g7_geometry_test_materials();
        n_azim=4, spacing=0.5, n_polar=2, bcs=bcs
    )
    sol = solve(prob; max_iterations=2, max_residual=0.0)

    @test length(tg.tracks_by_uid) > 0
    @test all(tag -> count(==(tag), prob.fsr_tag) > 0, 1:7)
    @test isfinite(sol.keff)
    @test sol.keff > 0
    @test isfinite(sol.residual)
    @test all(isfinite, sol.φ)
    @test sum(abs, sol.φ) > 0
end
