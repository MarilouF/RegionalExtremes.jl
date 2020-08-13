"""
    fit(m::Tuple{Integer, Integer}, p::Integer, y::Array{<:Real, 2},
        k::Integer, burnin::Integer, step::Vector{<:Real};
        μcov::Array{<:Real, 2} = zeros(size(g.y, 2), 0),
        ϕcov::Array{<:Real, 2} = zeros(size(g.y, 2), 0),
        ξcov::Array{<:Real, 2} = zeros(size(g.y, 2), 0),
        regional::Vector{Bool} = ones(Bool, 3 + size(μcov, 2) + size(ϕcov, 2) + size(ξcov, 2)))::Chain

Fit the generalized extreme value parameters for every grid cell using Gibbs sampling.

"""
function fit(m::Tuple{Integer, Integer}, p::Integer, y::Array{<:Real, 2},
    k::Integer, burnin::Integer, step::Vector{<:Real};
    μcov::Array{<:Real, 2} = zeros(size(y, 2), 0),
    ϕcov::Array{<:Real, 2} = zeros(size(y, 2), 0),
    ξcov::Array{<:Real, 2} = zeros(size(y, 2), 0),
    regional::Vector{Bool} = ones(Bool, 3 + size(μcov, 2) + size(ϕcov, 2) + size(ξcov, 2)))::Chain

    grid = BlockMaximaGrid(m, p, y)
    chain = Chain(grid, k, burnin, step, μcov = μcov, ϕcov = ϕcov, ξcov = ξcov, regional = regional)

    run!(chain)

    accrate = mean(chain.θacc[chain.burnin:end, :, :], dims=(1,3))[:, :, 1]'[:, 1]
    notrecommended = (0.2 .>= accrate) .| (accrate .>= 0.3)

    if any(notrecommended)
        @warn "Parameter(s) $(collect(1:size(chain.θacc, 2))[notrecommended]) acceptance rate $(round.(accrate[notrecommended], digits = 3)) not between 20% and 30%. Consider changing the step."
    end

    return chain
end
