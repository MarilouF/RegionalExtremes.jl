@testset "parameterestimation" begin
    include(joinpath("parameterestimation", "bayes_test.jl"))
    include(joinpath("parameterestimation", "fit_test.jl"))
end
