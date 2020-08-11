module RegionalExtremes
    using Distributions, ProgressMeter
    using Extremes, GMRF # TODO : Add GMRF as dependency

    include("util.jl")
    include("structures.jl")
    include("parameterestimation.jl")

    export BlockMaximaGrid, Chain, genGrid, fit
end # module
