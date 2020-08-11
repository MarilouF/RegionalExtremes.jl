@testset "fit.jl" begin
    @testset "fit(m, p, y, k, burnin, step; μcov, ϕcov, ξcov, regional)" begin
        # Calling fit calls run!
        chain = fit((2,2), 1, [1 0; 1 0; 1 0; 1 0], 4, 2, [0.5, 0.5, 0.5])

        # θ has values
        @test length(filter(p->p != 0, chain.θ)) > 0

        # θacc has values
        @test length(filter(p->p != 0, chain.θacc)) > 0

        # κ has values
        @test length(filter(p->p != 0, chain.κ)) > 0
    end
end
