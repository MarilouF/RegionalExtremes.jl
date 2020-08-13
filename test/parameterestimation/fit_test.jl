@testset "fit.jl" begin
    @testset "fit(m, p, y, k, burnin, step; μcov, ϕcov, ξcov, regional)" begin
        # Calling fit with inappropriate step warns
        chain = nothing
        @test_logs (:warn,) chain = RegionalExtremes.fit((2,2), 1, [1 0; 1 0; 1 0; 1 0], 4, 2, [0.1, 0.1, 0.1])

        # θ has values
        @test length(filter(p->p != 0, chain.θ)) > 0

        # θacc has values
        @test length(filter(p->p != 0, chain.θacc)) > 0

        # κ has values
        @test length(filter(p->p != 0, chain.κ)) > 0

        # Calling fit with appropriate step does not warn
        κs = [50.0, 50.0, 3000.0]
        g, θ, x = genGrid((10, 12), 1, κs..., 300)
        @test_logs RegionalExtremes.fit(g.G.gridSize, g.p, g.y, 500, 250, [0.2, 0.15, 0.08])
    end
end
