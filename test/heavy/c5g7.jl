# OECD/NEA C5G7 2D MOX benchmark MCNP reference values. The geometry and cross
# sections are specified in NEA/NSC/DOC(2001)4; the benchmark comparison report is
# NEA/NSC/DOC(2003)16 for deterministic transport without pin-cell homogenization.
const C5G7_2D_BENCHMARK_REFERENCE = (
    source="OECD/NEA C5G7 2D MOX benchmark MCNP reference",
    keff=1.18655,
    keff_uncertainty_pcm=9.5,
    max_normalized_pin_power=2.498,
    max_normalized_pin_power_uncertainty_percent=0.16,
    min_normalized_pin_power=0.232,
    min_normalized_pin_power_uncertainty_percent=0.58,
    assembly_powers=(
        lower_left_uo2=492.8,
        lower_right_mox=211.7,
        upper_left_mox=211.7,
        upper_right_uo2=139.8,
    ),
    assembly_power_uncertainty_percent=(
        lower_left_uo2=0.10,
        lower_right_mox=0.18,
        upper_left_mox=0.18,
        upper_right_uo2=0.20,
    ),
)

# OpenMOC does not commit a C5G7 `results_true.dat` in its regression suite.
# These values are the published OpenMOC C5G7 results from Boyd et al. (2014),
# as tabulated in Tramm et al. (2018), Table 2. The paper reports distribution
# error metrics against the official benchmark pin-power map, not absolute pin
# powers, so they are used as method-comparison guardrails below.
const C5G7_2D_OPENMOC_REFERENCE = (
    source="Boyd et al. (2014) OpenMOC C5G7, summarized by Tramm et al. (2018), Table 2",
    n_fsrs=142_964,
    n_azim=64,
    keff=1.18650,
    keff_error_pcm=-5.0,
    average_pin_power_error_percent=0.451,
    max_pin_power_error_percent=1.772,
    cmfd=(
        n_azim=64,
        keff=1.18659,
        keff_error_pcm=3.0,
        average_pin_power_error_percent=0.451,
        max_pin_power_error_percent=1.772,
    ),
)

const C5G7_2D_REFERENCE = C5G7_2D_BENCHMARK_REFERENCE

const C5G7_BENCHMARK_SETTINGS = (
    pin_lc=0.25,
    n_azim=32,
    spacing=0.1,
    n_polar=6,
    polar_quadrature=TabuchiYamamoto,
    segmentize_kwargs=(k=10, rtol=2 * sqrt(eps(Float64))),
    max_iterations=700,
    max_residual=1e-4,
    parallel=true,
    keff_tolerance_pcm=1500.0,
    material_volume_rtol=3e-2,
    pin_power_extrema_rtol=0.35,
    assembly_power_rtol=0.20,
    openmoc_keff_tolerance_pcm=1000.0,
    openmoc_pin_power_error_margin_percent=1.0,
)

const C5G7_HIGH_ACCURACY_BENCHMARK_OVERRIDES = (
    pin_lc=0.10,
    segmentize_kwargs=(k=20, rtol=2 * sqrt(eps(Float64))),
    keff_tolerance_pcm=150.0,
    material_volume_rtol=6e-3,
    pin_power_extrema_rtol=0.03,
    assembly_power_rtol=0.02,
    openmoc_keff_tolerance_pcm=150.0,
    openmoc_pin_power_error_margin_percent=0.25,
)

function c5g7_benchmark_settings(; kwargs...)
    return merge(C5G7_BENCHMARK_SETTINGS, NamedTuple(kwargs))
end

function c5g7_high_accuracy_benchmark_settings(; kwargs...)
    return merge(
        C5G7_BENCHMARK_SETTINGS,
        C5G7_HIGH_ACCURACY_BENCHMARK_OVERRIDES,
        NamedTuple(kwargs),
    )
end

function c5g7_2d_benchmark_result(; kwargs...)
    settings = c5g7_benchmark_settings(; kwargs...)
    reference = C5G7_2D_REFERENCE

    bcs = BoundaryConditions(top=Vacuum, bottom=Reflective, left=Reflective, right=Vacuum)
    geometry_time = @elapsed begin
        model = c5g7_benchmark_model(; pin_lc=settings.pin_lc)
    end
    track_generator_time = @elapsed begin
        tg = TrackGenerator(model, settings.n_azim, settings.spacing;
            bcs=bcs, volume_correction=true
        )
    end
    trace_time = @elapsed trace!(tg)
    segmentize_time = @elapsed segmentize!(tg; settings.segmentize_kwargs...)
    problem_time = @elapsed begin
        prob = MoCProblem(
            tg, settings.polar_quadrature(settings.n_polar), c5g7_benchmark_materials()
        )
    end
    solve_time = @elapsed begin
        sol = solve(prob;
            max_iterations=settings.max_iterations,
            max_residual=settings.max_residual,
            parallel=settings.parallel
        )
    end

    material_volumes = c5g7_material_volumes(model)
    reference_material_volumes = c5g7_analytic_material_volumes()
    power_result = c5g7_power_result(prob, sol)
    keff_error_pcm = 1.0e5 * (sol.keff - reference.keff)
    timings = (;
        geometry=geometry_time,
        track_generator=track_generator_time,
        trace=trace_time,
        segmentize=segmentize_time,
        problem=problem_time,
        solve=solve_time,
    )
    return (;
        settings,
        prob,
        tg,
        sol,
        power_result,
        keff_error_pcm,
        material_volumes,
        reference_material_volumes,
        timings,
        n_cells=length(tg.volumes),
        n_tracks=length(tg.tracks_by_uid),
        n_segments=count_track_segments(tg),
    )
end

function count_track_segments(tg)
    n_segments = 0
    for track in tg.tracks_by_uid
        n_segments += length(track.segments)
    end
    return n_segments
end

function c5g7_power_result(prob, sol)
    mesh = prob.trackgenerator.mesh
    NGroups = NeutronTransport.ngroups(prob)

    raw_pin_powers = zeros(4, C5G7_N_PINS, C5G7_N_PINS)
    raw_group_pin_powers = zeros(NGroups, 4, C5G7_N_PINS, C5G7_N_PINS)
    raw_assembly_powers = zeros(4)
    raw_group_assembly_powers = zeros(NGroups, 4)

    for cell in eachindex(prob.cell_to_fsr)
        fsr = prob.cell_to_fsr[cell]
        xs_tag = prob.fsr_tag[fsr]
        xs_tag in C5G7_ACTIVE_FUEL_TAGS || continue

        x, y = c5g7_cell_centroid(mesh, cell)
        pin_index = c5g7_pin_index(x, y)
        isnothing(pin_index) && throw(ArgumentError(
            "active C5G7 fuel cell $cell is outside the modeled fuel assemblies."
        ))
        assembly_idx, pin_i, pin_j = pin_index

        power = zero(eltype(prob))
        for g in 1:NGroups
            group_power = c5g7_fission_rate(prob, sol, fsr, g)
            raw_group_pin_powers[g, assembly_idx, pin_i, pin_j] += group_power
            raw_group_assembly_powers[g, assembly_idx] += group_power
            power += group_power
        end

        raw_pin_powers[assembly_idx, pin_i, pin_j] += power
        raw_assembly_powers[assembly_idx] += power
    end

    active_pin_mask = raw_pin_powers .> 0
    n_active_pins = count(active_pin_mask)
    mean_pin_power = sum(raw_pin_powers) / n_active_pins
    normalized_pin_powers = raw_pin_powers ./ mean_pin_power
    normalized_group_pin_powers = raw_group_pin_powers ./ mean_pin_power
    normalized_active_pin_powers = normalized_pin_powers[active_pin_mask]
    normalized_assembly_powers = raw_assembly_powers ./ mean_pin_power
    normalized_group_assembly_powers = raw_group_assembly_powers ./ mean_pin_power
    group_power_fractions = vec(sum(raw_group_pin_powers; dims=(2, 3, 4))) ./ sum(raw_pin_powers)

    assembly_powers = NamedTuple{C5G7_ASSEMBLY_LABELS}(normalized_assembly_powers)

    return (;
        raw_pin_powers,
        raw_group_pin_powers,
        normalized_pin_powers,
        normalized_group_pin_powers,
        active_pin_mask,
        normalized_active_pin_powers,
        assembly_powers,
        normalized_group_assembly_powers,
        group_power_fractions,
        max_normalized_pin_power=maximum(normalized_active_pin_powers),
        min_normalized_pin_power=minimum(normalized_active_pin_powers),
        n_active_pins,
    )
end

# Official C5G7 pin powers are fission-rate tallies. This follows OpenMOC's default
# `computeFSRFissionRates(..., nu=false)` convention: use Σf here, not νΣf.
function c5g7_fission_rate(prob, sol, fsr::Integer)
    NGroups = NeutronTransport.ngroups(prob)
    rate = zero(eltype(prob))
    @inbounds for g in 1:NGroups
        rate += c5g7_fission_rate(prob, sol, fsr, g)
    end
    return rate
end

function c5g7_fission_rate(prob, sol, fsr::Integer, g::Integer)
    xs = NeutronTransport.getxs(prob, fsr)
    volume = prob.volumes[fsr]
    NGroups = NeutronTransport.ngroups(prob)
    ig = (fsr - 1) * NGroups + g
    return xs.Σf[g] * sol.φ[ig] * volume
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
        Σf=[
            4.79002e-9, 5.82564e-9, 4.63719e-7, 5.24406e-6,
            1.45390e-7, 7.14972e-7, 2.08041e-6
        ],
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
        Σf=[
            0.00721206, 0.000819301, 0.00645320, 0.0185648,
            0.0178084, 0.0830348, 0.216004
        ],
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
        Σf=[
            0.00762704, 0.000876898, 0.00569835, 0.0228872,
            0.0107635, 0.232757, 0.248968
        ],
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
        Σf=[
            0.00825446, 0.00132565, 0.00842156, 0.0328730,
            0.0159636, 0.323794, 0.362803
        ],
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
        Σf=[
            0.00867209, 0.00162426, 0.0102716, 0.0390447,
            0.0192576, 0.374888, 0.430599
        ],
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

if !isdefined(@__MODULE__, :C5G7_RUN_TESTSET)
    C5G7_RUN_TESTSET = true
end

function c5g7_benchmark_errors(result)
    reference = C5G7_2D_REFERENCE
    power_result = result.power_result

    material_volume_relerr = abs.(
        (result.material_volumes .- result.reference_material_volumes) ./
        result.reference_material_volumes
    )
    reference_assembly_powers = reference.assembly_powers
    assembly_power_relerr = abs.(
        (collect(values(power_result.assembly_powers)) .-
         collect(values(reference_assembly_powers))) ./
        collect(values(reference_assembly_powers))
    )
    max_pin_power_relerr = abs(
        power_result.max_normalized_pin_power - reference.max_normalized_pin_power
    ) / reference.max_normalized_pin_power
    min_pin_power_relerr = abs(
        power_result.min_normalized_pin_power - reference.min_normalized_pin_power
    ) / reference.min_normalized_pin_power

    return (;
        material_volume_relerr,
        assembly_power_relerr,
        max_material_volume_relerr=maximum(material_volume_relerr),
        max_assembly_power_relerr=maximum(assembly_power_relerr),
        max_pin_power_relerr,
        min_pin_power_relerr,
        max_pin_power_extrema_relerr=max(max_pin_power_relerr, min_pin_power_relerr),
    )
end

function test_c5g7_benchmark_result(result)
    settings = result.settings
    reference = C5G7_2D_REFERENCE
    openmoc_reference = C5G7_2D_OPENMOC_REFERENCE
    prob = result.prob
    tg = result.tg
    sol = result.sol
    power_result = result.power_result
    errors = c5g7_benchmark_errors(result)
    openmoc_keff_error_pcm = 1.0e5 * (sol.keff - openmoc_reference.keff)
    openmoc_cmfd_keff_error_pcm = 1.0e5 * (sol.keff - openmoc_reference.cmfd.keff)
    max_pin_power_extrema_error_percent = 100 * errors.max_pin_power_extrema_relerr

    @info(
        "C5G7 2D benchmark comparison",
        benchmark_reference = reference.source,
        reference_keff = reference.keff,
        calculated_keff = sol.keff,
        keff_error_pcm = result.keff_error_pcm,
        reference_keff_uncertainty_pcm = reference.keff_uncertainty_pcm,
        max_material_volume_relerr = errors.max_material_volume_relerr,
        reference_max_pin_power = reference.max_normalized_pin_power,
        reference_max_pin_power_uncertainty_percent =
            reference.max_normalized_pin_power_uncertainty_percent,
        calculated_max_pin_power = power_result.max_normalized_pin_power,
        max_pin_power_relerr = errors.max_pin_power_relerr,
        reference_min_pin_power = reference.min_normalized_pin_power,
        reference_min_pin_power_uncertainty_percent =
            reference.min_normalized_pin_power_uncertainty_percent,
        calculated_min_pin_power = power_result.min_normalized_pin_power,
        min_pin_power_relerr = errors.min_pin_power_relerr,
        reference_assembly_powers = reference.assembly_powers,
        reference_assembly_power_uncertainty_percent =
            reference.assembly_power_uncertainty_percent,
        calculated_assembly_powers = power_result.assembly_powers,
        max_assembly_power_relerr = errors.max_assembly_power_relerr,
        openmoc_reference = openmoc_reference.source,
        openmoc_keff = openmoc_reference.keff,
        openmoc_keff_error_pcm,
        openmoc_cmfd_keff = openmoc_reference.cmfd.keff,
        openmoc_cmfd_keff_error_pcm,
        openmoc_average_pin_power_error_percent =
            openmoc_reference.average_pin_power_error_percent,
        openmoc_max_pin_power_error_percent = openmoc_reference.max_pin_power_error_percent,
        max_pin_power_extrema_error_percent,
        iterations = sol.iterations,
        residual = sol.residual,
        parallel = settings.parallel,
        threads = Base.Threads.nthreads(),
        pin_lc = settings.pin_lc,
        segmentize_kwargs = settings.segmentize_kwargs,
        n_cells = result.n_cells,
        n_tracks = result.n_tracks,
        n_segments = result.n_segments,
        timings = result.timings,
    )

    material_tags = Int32.(1:7)

    @test errors.max_material_volume_relerr <= settings.material_volume_rtol
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
    @test power_result.n_active_pins == C5G7_ACTIVE_FUEL_PIN_COUNT
    @test sum(power_result.normalized_active_pin_powers) ≈ C5G7_ACTIVE_FUEL_PIN_COUNT
    @test errors.max_pin_power_relerr <= settings.pin_power_extrema_rtol
    @test errors.min_pin_power_relerr <= settings.pin_power_extrema_rtol
    @test errors.max_assembly_power_relerr <= settings.assembly_power_rtol
    @test abs(openmoc_keff_error_pcm) <= settings.openmoc_keff_tolerance_pcm
    @test abs(openmoc_cmfd_keff_error_pcm) <= settings.openmoc_keff_tolerance_pcm

    # OpenMOC's published C5G7 table reports APE/MPE against the full MCNP pin-power
    # map. Until that full 34x34 pin map is encoded here, compare the strongest
    # official pin-power extrema error available in this test against OpenMOC's MPE.
    @test max_pin_power_extrema_error_percent <= (
        openmoc_reference.max_pin_power_error_percent +
        settings.openmoc_pin_power_error_margin_percent
    )

    return nothing
end

if C5G7_RUN_TESTSET
@testset "C5G7 fast 7-group solve" begin
    result = c5g7_2d_benchmark_result()
    test_c5g7_benchmark_result(result)
end

@testset "C5G7 high-accuracy 7-group solve" begin
    settings = c5g7_high_accuracy_benchmark_settings()
    result = c5g7_2d_benchmark_result(; settings...)
    test_c5g7_benchmark_result(result)
end
end
