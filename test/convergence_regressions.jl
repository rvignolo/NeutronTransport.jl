@testset "OpenMOC homogeneous discretization ladder" begin
    analytic_keff, analytic_flux_ratio = openmoc_homogeneous_medium_keff()
    material = openmoc_homogeneous_medium()
    model = openmoc_homogeneous_grid_model()

    configs = (
        (n_azim=4, spacing=0.45, n_polar=2),
        (n_azim=8, spacing=0.30, n_polar=2),
        (n_azim=8, spacing=0.20, n_polar=4),
    )

    keffs = Float64[]
    flux_ratios = Float64[]

    for config in configs
        prob = reflected_problem(model, [material];
            n_azim=config.n_azim,
            spacing=config.spacing,
            n_polar=config.n_polar
        )
        sol = solve(prob; max_iterations=300, max_residual=1e-8)
        group_1 = collect(sol(1))
        group_2 = collect(sol(2))

        push!(keffs, sol.keff)
        push!(flux_ratios, sum(group_2 ./ group_1) / length(group_1))
    end

    keff_errors = abs.(keffs .- analytic_keff)
    flux_ratio_errors = abs.(flux_ratios .- analytic_flux_ratio)

    @test all(<=(5e-4), keff_errors)
    @test all(<=(5e-4), flux_ratio_errors)
    @test last(keff_errors) <= 1.25 * first(keff_errors)
    @test last(flux_ratio_errors) <= 1.25 * first(flux_ratio_errors)
end
