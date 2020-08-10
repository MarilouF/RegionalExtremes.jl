module RegionalExtremes
    using Extremes, GMRF

    include("util.jl")
    include("structures.jl")
    include("parameterestimation.jl")

    export BlockMaximaGrid, Chain, genGrid, fit
end # module
