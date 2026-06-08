# OECD/NEA C5G7 2D MOX benchmark MCNP reference values. The current test asserts
# keff; pin and assembly powers are recorded here for later pin-power observables.
const C5G7_2D_REFERENCE = (
    keff=1.18655,
    keff_uncertainty=0.00006,
    max_normalized_pin_power=2.498,
    min_normalized_pin_power=0.232,
    assembly_powers=(492.8, 211.7, 139.8),
)

const C5G7_BENCHMARK_SETTINGS = (
    pin_lc=0.25,
    n_azim=32,
    spacing=0.1,
    n_polar=6,
    segmentize_kwargs=(k=10, rtol=2 * sqrt(eps(Float64))),
    max_iterations=700,
    max_residual=1e-4,
    keff_tolerance_pcm=1500.0,
    material_volume_rtol=3e-2,
)

function c5g7_2d_benchmark_result()
    settings = C5G7_BENCHMARK_SETTINGS
    reference = C5G7_2D_REFERENCE

    bcs = BoundaryConditions(top=Vacuum, bottom=Reflective, left=Reflective, right=Vacuum)
    model = c5g7_benchmark_model(; pin_lc=settings.pin_lc)
    tg = TrackGenerator(model, settings.n_azim, settings.spacing;
        bcs=bcs, volume_correction=true
    )
    trace!(tg)
    segmentize!(tg; settings.segmentize_kwargs...)
    prob = MoCProblem(tg, TabuchiYamamoto(settings.n_polar), c5g7_benchmark_materials())
    sol = solve(prob;
        max_iterations=settings.max_iterations,
        max_residual=settings.max_residual
    )

    material_volumes = c5g7_material_volumes(model)
    reference_material_volumes = c5g7_analytic_material_volumes()
    keff_error_pcm = 1.0e5 * (sol.keff - reference.keff)
    return (;
        prob, tg, sol, keff_error_pcm, material_volumes, reference_material_volumes
    )
end

function c5g7_benchmark_materials()
    n_groups = 7
    χ = [0.58791, 0.41176, 0.00033906, 1.1761e-7, 0.0, 0.0, 0.0]
    χ ./= sum(χ)

    guide_tube = CrossSections("guide-tube", n_groups;
        Σt=[0.126032, 0.29316, 0.28424, 0.28096, 0.33444, 0.56564, 1.17215],
        Σs0=[0.0661659 0.05907 0.00028334 1.4622e-6 2.0642e-8 0.0 0.0;
            0.0000000 0.240377 0.05243500 0.0002499 1.9239e-5 2.9875e-6 4.214e-7;
            0.0000000 0.000000 0.18329700 0.092397 0.0069446 0.0010803 0.00020567;
            0.0000000 0.000000 0.00000000 0.0788511 0.17014 0.025881 0.0049297;
            0.0000000 0.000000 0.00000000 3.7333e-5 0.0997372 0.20679 0.024478;
            0.0000000 0.000000 0.00000000 0.0 0.00091726 0.316765 0.23877;
            0.0000000 0.000000 0.00000000 0.0 0.0 0.049792 1.09912]
    )

    fission_chamber = CrossSections("fission-chamber", n_groups;
        χ=χ,
        νΣf=[1.323401e-8, 1.4345e-8, 1.128599e-6, 1.276299e-5, 3.538502e-7, 1.740099e-6, 5.063302e-6],
        Σt=[0.126032, 0.29316, 0.28425, 0.28102, 0.33446, 0.56564, 1.17214],
        Σs0=[0.0661659 0.05907 0.00028334 1.4622e-6 2.0642e-8 0.0 0.0;
            0.0 0.240377 0.052435 0.0002499 1.9239e-5 2.9875e-6 4.214e-7;
            0.0 0.0 0.183425 0.092288 0.0069365 0.001079 0.00020543;
            0.0 0.0 0.0 0.0790769 0.16999 0.02586 0.0049256;
            0.0 0.0 0.0 3.734e-5 0.099757 0.20679 0.024478;
            0.0 0.0 0.0 0.0 0.00091742 0.316774 0.23876;
            0.0 0.0 0.0 0.0 0.0 0.049793 1.0991]
    )

    uo2 = CrossSections("UO2", n_groups;
        χ=χ,
        νΣf=[0.02005998, 0.002027303, 0.01570599, 0.04518301, 0.04334208, 0.2020901, 0.5257105],
        Σt=[0.177949, 0.329805, 0.480388, 0.554367, 0.311801, 0.395168, 0.564406],
        Σs0=[0.127537 0.042378 9.4374e-6 5.5163e-9 0.0 0.0 0.0;
            0.0 0.324456 0.0016314 3.1427e-9 0.0 0.0 0.0;
            0.0 0.0 0.45094 0.0026792 0.0 0.0 0.0;
            0.0 0.0 0.0 0.452565 0.0055664 0.0 0.0;
            0.0 0.0 0.0 0.00012525 0.271401 0.010255 1.0021e-8;
            0.0 0.0 0.0 0.0 0.0012968 0.265802 0.016809;
            0.0 0.0 0.0 0.0 0.0 0.0085458 0.27308]
    )

    mox_43 = CrossSections("MOX_43", n_groups;
        χ=χ,
        νΣf=[0.021753, 0.002535103, 0.01626799, 0.0654741, 0.03072409, 0.666651, 0.7139904],
        Σt=[0.178731, 0.330849, 0.483772, 0.566922, 0.426227, 0.678997, 0.682852],
        Σs0=[0.128876 0.041413 8.229e-6 5.0405e-9 0.0 0.0 0.0;
            0.0 0.325452 0.0016395 1.5982e-9 0.0 0.0 0.0;
            0.0 0.0 0.453188 0.0026142 0.0 0.0 0.0;
            0.0 0.0 0.0 0.457173 0.0055394 0.0 0.0;
            0.0 0.0 0.0 0.00016046 0.276814 0.0093127 9.1656e-9;
            0.0 0.0 0.0 0.0 0.0020051 0.252962 0.01485;
            0.0 0.0 0.0 0.0 0.0 0.0084948 0.265007]
    )

    mox_7 = CrossSections("MOX_7", n_groups;
        χ=χ,
        νΣf=[0.02381395, 0.003858689, 0.024134, 0.09436622, 0.04576988, 0.9281814, 1.0432],
        Σt=[0.181323, 0.334368, 0.493785, 0.591216, 0.474198, 0.833601, 0.853603],
        Σs0=[0.130457 0.041792 8.5105e-6 5.1329e-9 0.0 0.0 0.0;
            0.0 0.328428 0.0016436 2.2017e-9 0.0 0.0 0.0;
            0.0 0.0 0.458371 0.0025331 0.0 0.0 0.0;
            0.0 0.0 0.0 0.463709 0.0054766 0.0 0.0;
            0.0 0.0 0.0 0.00017619 0.282313 0.0087289 9.0016e-9;
            0.0 0.0 0.0 0.0 0.002276 0.249751 0.013114;
            0.0 0.0 0.0 0.0 0.0 0.0088645 0.259529]
    )

    mox_87 = CrossSections("MOX_87", n_groups;
        χ=χ,
        νΣf=[0.025186, 0.004739509, 0.02947805, 0.11225, 0.05530301, 1.074999, 1.239298],
        Σt=[0.183045, 0.336705, 0.500507, 0.606174, 0.502754, 0.921028, 0.955231],
        Σs0=[0.131504 0.042046 8.6972e-6 5.1938e-9 0.0 0.0 0.0;
            0.0 0.330403 0.0016463 2.6006e-9 0.0 0.0 0.0;
            0.0 0.0 0.461792 0.0024749 0.0 0.0 0.0;
            0.0 0.0 0.0 0.468021 0.005433 0.0 0.0;
            0.0 0.0 0.0 0.00018597 0.285771 0.0083973 8.928e-9;
            0.0 0.0 0.0 0.0 0.0023916 0.247614 0.012322;
            0.0 0.0 0.0 0.0 0.0 0.0089681 0.256093]
    )

    water = CrossSections("water", n_groups;
        Σt=[0.159206, 0.41297, 0.59031, 0.58435, 0.718, 1.25445, 2.65038],
        Σs0=[0.0444777 0.1134 0.00072347 3.7499e-6 5.3184e-8 0.0 0.0;
            0.0 0.282334 0.12994 0.0006234 4.8002e-5 7.4486e-6 1.0455e-6;
            0.0 0.0 0.345256 0.22457 0.016999 0.0026443 0.00050344;
            0.0 0.0 0.0 0.0910284 0.41551 0.063732 0.012139;
            0.0 0.0 0.0 7.1437e-5 0.139138 0.51182 0.061229;
            0.0 0.0 0.0 0.0 0.0022157 0.699913 0.53732;
            0.0 0.0 0.0 0.0 0.0 0.13244 2.4807]
    )

    return [guide_tube, fission_chamber, uo2, mox_43, mox_7, mox_87, water]
end

@testset "C5G7 full 7-group demo solve" begin
    settings = C5G7_BENCHMARK_SETTINGS
    reference = C5G7_2D_REFERENCE
    result = c5g7_2d_benchmark_result()
    prob = result.prob
    tg = result.tg
    sol = result.sol

    material_volume_relerr = abs.(
        (result.material_volumes .- result.reference_material_volumes) ./
        result.reference_material_volumes
    )

    @info(
        "C5G7 2D benchmark comparison",
        reference_keff = reference.keff,
        calculated_keff = sol.keff,
        keff_error_pcm = result.keff_error_pcm,
        max_material_volume_relerr = maximum(material_volume_relerr),
    )

    material_tags = Int32.(1:7)

    @test maximum(material_volume_relerr) <= settings.material_volume_rtol
    @test length(tg.tracks_by_uid) > 0
    @test sort(unique(prob.fsr_tag)) == material_tags
    @test length(prob.fsr_tag) == length(prob.volumes)
    @test length(prob.cell_to_fsr) == length(tg.volumes)
    @test all(>(0), prob.volumes)
    @test sum(prob.volumes) ≈ sum(tg.volumes)
    @test all(tag -> sum(prob.volumes[prob.fsr_tag .== tag]) > 0, material_tags)
    @test NeutronTransport.ngroups(prob) == 7
    @test isfinite(sol.keff)
    @test isfinite(sol.residual)
    @test all(isfinite, sol.φ)
    @test sum(abs, sol.φ) > 0
    @test sol.iterations < settings.max_iterations
    @test sol.residual <= settings.max_residual
    @test abs(result.keff_error_pcm) <= settings.keff_tolerance_pcm
end
