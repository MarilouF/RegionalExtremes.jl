using RegionalExtremes
using Test

@testset "RegionalExtremes.jl" begin
    include("util_test.jl")
    include("structures_test.jl")
    include("parameterestimation_test.jl")
end
