@testset "polar quadrature constructors" begin
    for constructor in (
        TabuchiYamamoto,
        GaussLegendre,
        EqualWeight,
        EqualAngle,
        Leonard,
    )
        for n_polar in (2, 4, 6)
            quad = constructor(n_polar)
            half = div(n_polar, 2)
            lower = 1:half
            upper = half+1:n_polar

            @test NeutronTransport.n_polar(quad) == n_polar
            @test NeutronTransport.n_polar_half(quad) == half
            @test all(isfinite, quad.sinθs)
            @test all(isfinite, quad.θs)
            @test all(isfinite, quad.ωₚ)
            @test quad.sinθs ≈ sin.(quad.θs)
            @test sum(quad.ωₚ) ≈ 1
            @test quad.sinθs[lower] ≈ reverse(quad.sinθs[upper])
            @test quad.ωₚ[lower] ≈ reverse(quad.ωₚ[upper])
        end
    end
end
