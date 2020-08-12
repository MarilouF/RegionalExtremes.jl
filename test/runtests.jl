using RegionalExtremes
using Test, Statistics, Distributions
import Random

Random.seed!(2)

@testset "RegionalExtremes.jl" begin
    include("util_test.jl")
    include("structures_test.jl")
    include("parameterestimation_test.jl")
end
