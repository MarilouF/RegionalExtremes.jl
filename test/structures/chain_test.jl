@testset "chain.jl" begin
    @testset "Chain(g, k, burnin, step; μcov, ϕcov, ξcov, regional)" begin
        # TODO : Test size(μcov, 1) != n
        # TODO : Test size(ϕcov, 1) != n
        # TODO : Test size(ξcov, 1) != n
        # TODO : length(step) != np
        # TODO : length(regional) != np
        # TODO : Test all Chain attributes
    end
    @testset "gev_loglikelihood(chain, candidate, a)" begin
        # TODO : Test stationary with known values
        # TODO : Test non stationary with known values
    end
    @testset "local_loglikelihood(chain, candidates, subset, i)" begin
        # TODO : Test with known values subset size 1
    end
    @testset "regional_loglikelihood(chain, candidates, subset, i, param)" begin
        # TODO : Test with known values subset size 1
    end
end
