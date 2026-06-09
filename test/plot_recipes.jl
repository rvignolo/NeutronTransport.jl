@testset "plot recipe data helpers" begin
    materials = [openmoc_homogeneous_medium(), openmoc_water_medium()]
    bcs = BoundaryConditions(
        top=Reflective, bottom=Reflective, left=Reflective, right=Reflective
    )
    tg = TrackGenerator(openmoc_reflective_grid_model(), 8, 0.25;
        bcs=bcs, volume_correction=true
    )
    trace!(tg)
    segmentize!(tg)

    cell_to_fsr = cell_material_ids(tg, materials)
    prob = MoCProblem(tg, TabuchiYamamoto(2), materials; cell_to_fsr=cell_to_fsr)
    sol = solve(prob; max_iterations=2, max_residual=0.0)

    cell_flux = cell_scalar_flux(sol, 1)
    @test length(cell_flux) == length(prob.cell_to_fsr)
    @test cell_flux ≈ fsr_to_cell_values(prob, sol(1))

    for cell in eachindex(cell_flux)
        @test cell_flux[cell] == sol(Int(prob.cell_to_fsr[cell]), 1)
    end

    field = CellScalarField(sol, 1)
    @test field.mesh === tg.mesh
    @test field.values == cell_flux
    @test field.title == "Scalar flux group 1"

    @test_throws ArgumentError fsr_to_cell_values(prob, sol.φ)
    @test_throws ArgumentError CellScalarField(tg.mesh, cell_flux[begin:end-1])

    values = [1.0 2.0; 3.0 4.0]
    active = Bool[1 0; 0 1]
    pin_map = PinPowerMap(values; active=active, title=:Power)
    @test pin_map.values === values
    @test pin_map.active === active
    @test pin_map.title == "Power"
    @test_throws ArgumentError PinPowerMap(values; active=trues(1, 2))
end
