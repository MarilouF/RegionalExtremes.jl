@testset "bayes.jl" begin
    @testset "setinitialvalues!(chain)" begin
        # TODO : Test gumbel initial values with known values
        # TODO : Test κ initial values with known values (local and regional)
        # TODO : Test logpdf initial values with known values (local and regional)
    end
    @testset "copylast(chain, i)" begin
        # TODO : Test θ
        # TODO : Test κ
        # TODO : Test logpdf
    end
    @testset "candidates(chain, i, subset, param)" begin
        # TODO : Test only param column was changed
        # TODO : Test insupport
    end
    @testset "MHupdatechain!(chain, i, subset, param, c, logpdfP)" begin
        # TODO : Test all updated
        # TODO : Test all kept old
    end
    @testset "MHlocal!(chain, i, p)" begin
        # TODO : Check all iterations except i are not changed
    end
    @testset "MHregional!(chain, i, p)" begin
        # TODO : Check all iterations except i are not changed
    end
    @testset "MHrealisation!(chain, i)" begin
        # TODO : Check all iterations except i are not changed (local and regional)
    end
    @testset "κrealisation!(chain, i)" begin
        # TODO : Test in support with known values
    end
    @testset "run!(chain)" begin
        # TODO : Test with known values... (long tests will be here)
    end
end
