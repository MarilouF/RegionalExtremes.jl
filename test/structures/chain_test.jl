@testset "chain.jl" begin
    @testset "Chain(g, k, burnin, step; μcov, ϕcov, ξcov, regional)" begin
        m = (1, 2)
        p = 1
        y = [1 0; 1 0]
        k = 3000
        burnin = 500
        step = [0.5, 0.5, 0.5]

        # μcov not same length as data throws
        @test_throws AssertionError Chain(BlockMaximaGrid(m, p, y), k, burnin, step, μcov = zeros(0, 0))

        # ϕcov not same length as data throws
        @test_throws AssertionError Chain(BlockMaximaGrid(m, p, y), k, burnin, step, ϕcov = zeros(0, 0))

        # ξcov not same length as data throws
        @test_throws AssertionError Chain(BlockMaximaGrid(m, p, y), k, burnin, step, ξcov = zeros(0, 0))

        # step not same length as the number of parameters throws
        @test_throws AssertionError Chain(BlockMaximaGrid(m, p, y), k, burnin, [0])

        # regional not same length as the number of parameters throws
        @test_throws AssertionError Chain(BlockMaximaGrid(m, p, y), k, burnin, step, regional = zeros(Bool, 0))

        # Chain attributes
        step = [0.5, 0.5, 0.5, 0.5]
        regional = [true, true, true, false]
        μcov = ones(2, 1)
        chain = Chain(BlockMaximaGrid(m, p, y), k, burnin, step, μcov = μcov, regional = regional)

        @test chain.g.p == p
        @test chain.g.G.gridSize == m
        @test chain.k == k
        @test chain.burnin == burnin
        @test size(chain.θ) == (k, 4, prod(m))
        @test size(chain.θacc) == (k, 4, prod(m))
        @test chain.θstep == step
        @test size(chain.θlogpdf) == (k, 4, prod(m))
        @test chain.θisRegional == regional
        @test size(chain.κ) == (k, 3)
        @test chain.cov == [hcat(ones(2 ,1), μcov), ones(2, 1), ones(2, 1)]
        @test chain.index == [[1, 2], [3], [4]]
    end
    @testset "Base.show(io, obj)" begin
        io = IOBuffer()
        g = BlockMaximaGrid((1, 2), 1, [1 2 3; 4 5 6])
        c = Chain(g, 3000, 500, [0.5, 0.2, 0.1])
        @test_logs Base.show(io, c)
    end
    @testset "gev_loglikelihood(chain, candidate, a)" begin
        # Stationary
        chain = Chain(BlockMaximaGrid((1, 2), 1, [1 1; 1 1]), 3000, 500, [0.5, 0.5, 0.5])
        @test RegionalExtremes.gev_loglikelihood(chain, [1, 0, 0], 1) ≈ -2.0 atol = 0.001

        # Non-stationary
        μcov = zeros(2, 1)
        μcov[2, 1] = 1.0
        chain = Chain(BlockMaximaGrid((1, 2), 1, [1 2; 1 2]), 3000, 500, [0.5, 0.5, 0.5, 0.5], μcov = μcov)
        @test RegionalExtremes.gev_loglikelihood(chain, [1, 1, 0, 0], 1) ≈ -2.0 atol = 0.001
    end
    @testset "local_loglikelihood(chain, candidates, subset)" begin
        # Stationary with subset of size 1
        chain = Chain(BlockMaximaGrid((1, 2), 1, [1 1; 1 1]), 3000, 500, [0.5, 0.5, 0.5])
        @test RegionalExtremes.local_loglikelihood(chain, [1 1; 0 0; 0 0], [1])[] ≈ -2.0 atol = 0.001
    end
    @testset "regional_loglikelihood(chain, candidates, subset, i, param)" begin
        # Stationary with subset of size 1
        chain = Chain(BlockMaximaGrid((1, 2), 1, [1 1; 1 1]), 3000, 500, [0.5, 0.5, 0.5])
        chain.κ[1, 1] = 1.0
        @test RegionalExtremes.regional_loglikelihood(chain, [1 1; 0 0; 0 0], [1], 1, 1)[] ≈ -0.9189 atol = 0.001
    end
end
