function fit(m::Tuple{Integer, Integer}, p::Integer, y::Array{<:Real, 2},
        k::Integer, burnin::Integer, step::Vector{<:Real};
        μcov::Array{<:Real, 2} = zeros(size(g.y, 2), 0),
        ϕcov::Array{<:Real, 2} = zeros(size(g.y, 2), 0),
        ξcov::Array{<:Real, 2} = zeros(size(g.y, 2), 0),
        regional::Vector{Bool} = ones(Bool, 3 + size(μcov, 2) + size(ϕcov, 2) + size(ξcov, 2)))::Chain

        grid = BlockMaximaGrid(m, p, y)
        chain = Chain(grid, k, burnin, step, μcov = μcov, ϕcov = ϕcov, ξcov = ξcov, regional = regional)

        fitbayes!(chain)

        return chain
end
