@testset "CrossSections constructor" begin
    mixed = CrossSections("mixed", 2;
        Σt=Float64[0.4, 0.7],
        Σs0=Float32[0.3 0.05; 0.0 0.4]
    )

    @test eltype(mixed) === Float64
    @test mixed.Σt isa Vector{Float64}
    @test mixed.Σs0 isa Matrix{Float64}
    @test isapprox(mixed.Σs0_sum, [0.35, 0.4]; rtol=1e-6)
    @test mixed.χ == zeros(2)
    @test !NeutronTransport.isfissionable(mixed)

    fuel = CrossSections("fuel", 2;
        νΣf=[0.02, 0.3],
        Σt=[0.4, 0.7],
        Σs0=[0.3 0.05; 0.0 0.4]
    )

    @test NeutronTransport.isfissionable(fuel)
    @test fuel.χ == [1.0, 0.0]

    normalized = CrossSections("normalized", 2;
        χ=[0.999999999, 1e-9],
        νΣf=[0.02, 0.3],
        Σt=[0.4, 0.7],
        Σs0=[0.3 0.05; 0.0 0.4]
    )

    @test sum(normalized.χ) ≈ 1

    Σt_static = SVector(0.4, 0.7)
    Σs0_static = @SMatrix [0.3 0.05; 0.0 0.4]
    static = CrossSections("static", 2;
        νΣf=SVector(0.02, 0.3),
        Σt=Σt_static,
        Σs0=Σs0_static
    )

    @test static.χ isa SVector{2,Float64}
    @test static.Σt isa SVector{2,Float64}
    @test static.Σs0 isa SMatrix{2,2,Float64}
    @test static.Σs0_sum isa SVector{2,Float64}
    @test isapprox(static.Σs0_sum, SVector(0.35, 0.4); rtol=1e-6)

    @test_throws ArgumentError CrossSections("bad", 0; Σt=[1.0], Σs0=reshape([1.0], 1, 1))
    @test_throws ArgumentError CrossSections("bad", 2; Σt=[1.0], Σs0=[1.0 0.0; 0.0 1.0])
    @test_throws ArgumentError CrossSections("bad", 2; Σt=[1.0, 1.0], Σs0=[1.0 0.0])
    @test_throws ArgumentError CrossSections("bad", 2;
        χ=[0.4, 0.4],
        νΣf=[0.02, 0.3],
        Σt=[0.4, 0.7],
        Σs0=[0.3 0.05; 0.0 0.4]
    )
end
