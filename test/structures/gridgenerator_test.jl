@testset "gridgenerator.jl" begin
    @testset "genVar(t, n)" begin
        n = 1000

        # Generate rdm variable
        v = RegionalExtremes.genVar(RegionalExtremes.rdm, n)
        @test all(insupport.(Uniform(), v))

        # Generate temporal variable
        v = RegionalExtremes.genVar(RegionalExtremes.time, n)
        @test v == collect(1:n)
    end
    @testset "genParam(m, p, κ, n, ev, δ)" begin
        m = (1, 2)
        p = 1
        κ = 1.0
        n = 1000

        # Stationary
        ev = Vector{Tuple{RegionalExtremes.VarType, Real}}()
        np = length(ev)
        x, β₀, βᵢs, θ = RegionalExtremes.genParam(m, p, κ, n, ev)
        @test size(x) == (n, np)
        @test length(β₀) == prod(m)
        @test mean(β₀) ≈ 0 atol = 0.001
        @test size(βᵢs) == (prod(m), np)
        @test size(θ) == (prod(m), n)

        # Non-stationary
        ev = [(RegionalExtremes.rdm, 1.0)]
        np = length(ev)
        x, β₀, βᵢs, θ = RegionalExtremes.genParam(m, p, κ, n, ev)
        @test size(x) == (n, np)
        @test length(β₀) == prod(m)
        @test mean(β₀) ≈ 0.0 atol = 0.001
        @test size(βᵢs) == (prod(m), np)
        @test size(θ) == (prod(m), n)

        # δ is mean of generated parameter
        δ = 10
        x, β₀, βᵢs, θ = RegionalExtremes.genParam(m, p, κ, n, ev, δ)
        @test mean(β₀) ≈ 10.0 atol = 0.001
    end
    @testset "genGrid(m, p, κμ, κϕ, κξ, n; evμ, evϕ, evξ)" begin
        m = (1, 2)
        p = 1
        κs = [50.0, 50.0, 3000.0]
        n = 300

        # Stationary
        g, θ, x = genGrid(m, p, κs..., n)

        @test g.p == p
        @test g.G.gridSize == m
        @test size(g.y) == (prod(m), n)
        @test size(θ) == (prod(m), 3)
        @test size(x[1]) == (n, 0)
        @test size(x[2]) == (n, 0)
        @test size(x[3]) == (n, 0)

        # Non-stationary
        evμ = [(RegionalExtremes.rdm, 1.0), (RegionalExtremes.time, 1.0)]
        evϕ = [(RegionalExtremes.rdm, 1.0)]
        g, θ, x = genGrid(m, p, κs..., n, evμ = evμ, evϕ = evϕ)

        @test g.p == p
        @test g.G.gridSize == m
        @test size(g.y) == (prod(m), n)
        @test size(θ) == (prod(m), 3 + length(evμ) + length(evϕ))
        @test size(x[1]) == (n, length(evμ))
        @test size(x[2]) == (n, length(evϕ))
        @test size(x[3]) == (n, 0)
    end
end
