@testset "OpenMOC homogeneous infinite medium golden reference" begin
    prob, _ = single_cell_problem(BoundaryConditions(
        top=Reflective, bottom=Reflective, left=Reflective, right=Reflective
    ))
    reference = openmoc_homogeneous_reference()
    analytic_keff, analytic_flux_ratio = openmoc_homogeneous_medium_keff()
    prob = MoCProblem(prob.trackgenerator, TabuchiYamamoto(2), [openmoc_homogeneous_medium()])
    sol = solve(prob; max_iterations=300, max_residual=1e-8)

    @test sol.keff ≈ reference.keff rtol = 5e-4
    @test sol.keff ≈ analytic_keff rtol = 5e-4
    @test sol(1, 2) / sol(1, 1) ≈ reference.flux_ratio rtol = 5e-4
    @test sol(1, 2) / sol(1, 1) ≈ analytic_flux_ratio rtol = 5e-4
end

@testset "OpenMOC homogeneous 10x10 grid golden reference" begin
    reference = openmoc_homogeneous_reference()
    prob = reflected_problem(openmoc_homogeneous_grid_model(), [openmoc_homogeneous_medium()])
    sol = solve(prob; max_iterations=300, max_residual=1e-8)

    group_1 = collect(sol(1))
    group_2 = collect(sol(2))

    @test sol.keff ≈ reference.keff rtol = 5e-4
    flux_ratios = group_2 ./ group_1
    @test sum(flux_ratios) / length(flux_ratios) ≈ reference.flux_ratio rtol = 5e-4
    @test maximum(group_1) - minimum(group_1) <= 1e-7 * maximum(group_1)
    @test maximum(group_2) - minimum(group_2) <= 1e-7 * maximum(group_2)
end

@testset "OpenMOC pin-cell golden reference" begin
    reference = openmoc_pin_cell_reference()
    prob, tg = openmoc_pin_cell_problem(; coalesce_materials=true)
    sol = solve(prob; max_iterations=500, max_residual=1e-5)

    @test prob.fsr_tag == Int32[1, 2]
    @test count(==(1), prob.cell_to_fsr) > 0
    @test count(==(2), prob.cell_to_fsr) > 0
    @test sum(prob.volumes) ≈ 16.0
    @test sum(tg.volumes) ≈ 16.0
    @test isfinite(sol.keff)
    @test isfinite(sol.residual)
    @test all(isfinite, sol.φ)

    # The coalesced two-region FSR map matches OpenMOC's pin-cell setup.
    @test sol.keff ≈ reference.keff rtol = 5e-3
end

@testset "OpenMOC-style reflected 3x3 grid structural parity" begin
    prob = reflected_problem(
        openmoc_reflective_grid_model(),
        [openmoc_homogeneous_medium(), openmoc_water_medium()]
    )
    sol = solve(prob; max_iterations=20, max_residual=0.0)

    @test count(==(1), prob.fsr_tag) == 1
    @test count(==(2), prob.fsr_tag) == 8
    @test isfinite(sol.keff)
    @test isfinite(sol.residual)
    @test all(isfinite, sol.φ)
end

@testset "OpenMOC-style reflected 3x3 lattice-equivalent structural parity" begin
    materials = [openmoc_homogeneous_medium(), openmoc_water_medium()]
    grid_prob = reflected_problem(openmoc_reflective_grid_model(), materials)
    lattice_prob = reflected_problem(openmoc_lattice_grid_model(), materials)

    grid_sol = solve(grid_prob; max_iterations=20, max_residual=0.0)
    lattice_sol = solve(lattice_prob; max_iterations=20, max_residual=0.0)

    @test grid_prob.fsr_tag == lattice_prob.fsr_tag
    @test lattice_sol.keff ≈ grid_sol.keff
    @test lattice_sol.φ ≈ grid_sol.φ
end

@testset "OpenMOC-style reflected 3x2 oblong lattice structural parity" begin
    prob = reflected_problem(
        openmoc_oblong_lattice_grid_model(),
        [
            openmoc_homogeneous_medium(),
            openmoc_water_medium(),
            openmoc_guide_tube_medium(),
        ]
    )
    sol = solve(prob; max_iterations=20, max_residual=0.0)

    @test count(==(1), prob.fsr_tag) == 2
    @test count(==(2), prob.fsr_tag) == 2
    @test count(==(3), prob.fsr_tag) == 2
    @test isfinite(sol.keff)
    @test isfinite(sol.residual)
    @test all(isfinite, sol.φ)
end

@testset "OpenMOC-style 4x4 simple lattice material parity" begin
    expected = openmoc_simple_lattice_expected_volumes()
    prob, tg = openmoc_simple_lattice_problem(; coalesce_materials=true)
    sol = solve(prob; max_iterations=120, max_residual=1e-5)

    fuel_fsr = findfirst(==(Int32(1)), prob.fsr_tag)
    water_fsr = findfirst(==(Int32(2)), prob.fsr_tag)

    @test prob.fsr_tag == Int32[1, 2]
    @test count(==(1), prob.cell_to_fsr) > 0
    @test count(==(2), prob.cell_to_fsr) > 0
    @test prob.volumes[fuel_fsr] ≈ expected.fuel rtol = 5e-2
    @test prob.volumes[water_fsr] ≈ expected.water rtol = 5e-2
    @test sum(prob.volumes) ≈ expected.total
    @test sum(tg.volumes) ≈ expected.total
    @test isfinite(sol.keff)
    @test isfinite(sol.residual)
    @test all(isfinite, sol.φ)
end
