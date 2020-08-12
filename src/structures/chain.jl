struct Chain # TODO : println
    g::BlockMaximaGrid
    k::Integer
    burnin::Integer
    θ::Array{<:Real, 3}
    θacc::Array{Bool, 3}
    θstep::Vector{Real}
    θlogpdf::Array{<:Real, 3}
    θisRegional::Vector{Bool}
    κ::Array{Real, 2}
    cov::Vector{Array{T, 2}} where T<:Real
    index::Vector{Vector{T}} where T<:Integer
end

Base.Broadcast.broadcastable(obj::Chain) = Ref(obj)

"""
    Chain(g::BlockMaximaGrid, k::Integer, burnin::Integer,
        step::Vector{<:Real};
        μcov::Array{<:Real, 2} = zeros(size(g.y, 2), 0),
        ϕcov::Array{<:Real, 2} = zeros(size(g.y, 2), 0),
        ξcov::Array{<:Real, 2} = zeros(size(g.y, 2), 0),
        regional::Vector{Bool} = ones(Bool, 3 + size(μcov, 2) + size(ϕcov, 2) + size(ξcov, 2)))::Chain

Build a MCMC Chain.

"""
function Chain(g::BlockMaximaGrid, k::Integer, burnin::Integer,
    step::Vector{<:Real};
    μcov::Array{<:Real, 2} = zeros(size(g.y, 2), 0),
    ϕcov::Array{<:Real, 2} = zeros(size(g.y, 2), 0),
    ξcov::Array{<:Real, 2} = zeros(size(g.y, 2), 0),
    regional::Vector{Bool} = ones(Bool, 3 + size(μcov, 2) + size(ϕcov, 2) + size(ξcov, 2)))::Chain

    n = size(g.y, 2)
    @assert size(μcov, 1) == n "The number of covariate data for μ should correspond to the number of data per grid cell"
    @assert size(ϕcov, 1) == n "The number of covariate data for ϕ should correspond to the number of data per grid cell"
    @assert size(ξcov, 1) == n "The number of covariate data for ξ should correspond to the number of data per grid cell"

    np = 3 + size(μcov, 2) + size(ϕcov, 2) + size(ξcov, 2)
    @assert length(step) == np "The size of step should correspond to the number of parameters to estimate"
    @assert length(regional) == np "The size of regional should correspond to the number of parameters to estimate"

    m = prod(g.G.gridSize)

    θ = zeros(Float64, k, np, m)
    θacc = zeros(Int64, k, np, m)
    θlogpdf = zeros(Float64, k, np, m)
    κ = zeros(Float64, k, sum(regional))

    identitycol = ones(n, 1)
    cov = [hcat(identitycol, μcov), hcat(identitycol, ϕcov), hcat(identitycol, ξcov)]

    i = 0
    function increasei() # TODO : Find other solution
        i = i + 1
        return i
    end
    index = [
        Integer[increasei() for k in 1:size(μcov, 2) + 1],
        Integer[increasei() for k in 1:size(ϕcov, 2) + 1],
        Integer[increasei() for k in 1:size(ξcov, 2) + 1]
    ]

    return Chain(g, k, burnin, θ, θacc, step, θlogpdf, regional, κ, cov, index)
end

"""
    gev_loglikelihood(chain::Chain, candidate::Vector{<:Real}, a::Integer)::Real

Calculate the log likelihood for the candidates at grid cell `a` taking into account the covariates.

"""
function gev_loglikelihood(chain::Chain, candidate::Vector{<:Real}, a::Integer)::Real
    μ = multiply(chain.cov[1], candidate[chain.index[1]])
    ϕ = multiply(chain.cov[2], candidate[chain.index[2]])
    ξ = multiply(chain.cov[3], candidate[chain.index[3]])

    return sum(logpdf.(GeneralizedExtremeValue.(μ, exp.(ϕ), ξ), chain.g.y[a, :]))
end

"""
    local_loglikelihood(chain::Chain, candidates::Array{<:Real, 2}, subset::Vector{<:Integer}, i::Integer)::Vector{<:Real}

Calculate the log likelihood for the candidates at each individual grid cells in the subset.

"""
function local_loglikelihood(chain::Chain, candidates::Array{<:Real, 2}, subset::Vector{<:Integer}, i::Integer)::Vector{<:Real} # TODO : Remove i
    return gev_loglikelihood.(chain, [candidates[:, subset[i]] for i in 1:length(subset)], subset)
end

"""
    regional_loglikelihood(chain::Chain, candidates::Array{<:Real, 2}, subset::Vector{<:Integer}, i::Integer, param::Integer)::Vector{<:Real}

Calculate the Gaussian Marcov random fields full conditional log likelihood for the candidates for the grid cell in `subset` for the parameter `param`.

"""
function regional_loglikelihood(chain::Chain, candidates::Array{<:Real, 2}, subset::Vector{<:Integer}, i::Integer, param::Integer)::Vector{<:Real}
    return GMRF.fullcondlogpdf(GMRF.iGMRF(chain.g.G.gridSize..., chain.g.p, chain.κ[i, sum(chain.θisRegional[1:param])]), candidates[param, :])[subset]
end
