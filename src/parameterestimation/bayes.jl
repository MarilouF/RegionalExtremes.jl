"""
    setinitialvalues!(chain::Chain)

Set the initial values for the generalized extreme value and precision parameters.

"""
function setinitialvalues!(chain::Chain)
    gumbels = Extremes.gumbelfitpwm.([chain.g.y[i, :] for i in 1:size(chain.g.y, 1)])
    for i in 1:length(gumbels)
        for p in 1:3
            chain.θ[1, chain.index[p][1], i] = gumbels[i].θ̂[p]
        end
    end
    for p in 1:size(chain.θ, 2)
        chain.θlogpdf[1, p, :] = local_loglikelihood(chain, chain.θ[1, :, :], collect(1:prod(chain.g.G.gridSize)))
        if chain.θisRegional[p]
            κ = chain.θ[1, p, :]' * chain.g.G.W * chain.θ[1, p, :] ./ (prod(chain.g.G.gridSize) - 1)
            chain.κ[1, sum(chain.θisRegional[1:p])] = κ == 0 ? 1.0 : κ

            chain.θlogpdf[1, p, :] .+= regional_loglikelihood(chain, chain.θ[1, :, :], collect(1:prod(chain.g.G.gridSize)), 1, p)
        end
    end
end

"""
    copylast!(chain::Chain, i::Integer)

Copy the last Gibbs iteration into current iteration.

"""
function copylast!(chain::Chain, i::Integer)
    chain.θ[i, :, :] = chain.θ[i - 1, :, :]
    chain.κ[i, :] = chain.κ[i - 1, :]
    chain.θlogpdf[i, :, :] = chain.θlogpdf[i - 1, :, :]
end

"""
    candidates(chain::Chain, i::Integer, subset::Vector{<:Integer}, param::Integer)::Array{<:Real, 2}

Generate random step candidates for the parameter `param` and the subset of grid cell `subset`.

"""
function candidates(chain::Chain, i::Integer, subset::Vector{<:Integer}, param::Integer)::Array{<:Real, 2}
    c = copy(chain.θ[i, :, :])
    c[param, subset] = rand.(Normal.(chain.θ[i, param, subset], chain.θstep[param]))
    return c
end

"""
    MHupdatechain!(chain::Chain, i::Integer, subset::Vector{<:Integer},
        param::Integer, c::Array{<:Real, 2}, logpdfP::Vector{<:Real})

Update the chain for the accepted Metropolis Hastings candidates.

"""
function MHupdatechain!(chain::Chain, i::Integer, subset::Vector{<:Integer},
    param::Integer, c::Array{<:Real, 2}, logpdfP::Vector{<:Real})

    logα = min.(0, logpdfP - chain.θlogpdf[i, param, subset])
    logu = log(rand(Uniform()))

    a = logu .< logα

    chain.θacc[i, param, subset] = a
    chain.θ[i, param, subset[a]] = c[param, subset[a]]
    chain.θlogpdf[i, param, subset[a]] = logpdfP[a]
end

"""
    MHlocal!(chain::Chain, i::Integer, p::Integer)

Realize a Metropolis Hastings step for the parameter `p` for every grid cell without using a regional prior.

"""
function MHlocal!(chain::Chain, i::Integer, p::Integer)
    set = collect(1:size(chain.g.y, 1))
    c = candidates(chain, i, set, p)
    logpdf = local_loglikelihood(chain, c, set)
    MHupdatechain!(chain, i, set, p, c, logpdf)
end

"""
    MHregional!(chain::Chain, i::Integer, p::Integer)

Realize a Metropolis Hastings step for the parameter `p` for every grid cell using a regional prior.

"""
function MHregional!(chain::Chain, i::Integer, p::Integer)
    for subset in chain.g.G.condIndSubset
        c = candidates(chain, i, subset, p)
        logpdf = local_loglikelihood(chain, c, subset) .+ regional_loglikelihood(chain, c, subset, i, p)
        MHupdatechain!(chain, i, subset, p, c, logpdf)
    end
end

"""
    MHrealisation!(chain::Chain, i::Integer)

Realize a Metropolis Hastings step for every generalized extreme value parameter for every grid cell.

"""
function MHrealisation!(chain::Chain, i::Integer)
    for p in 1:size(chain.θ, 2)
        if chain.θisRegional[p]
            MHregional!(chain, i, p)
        else
            MHlocal!(chain, i, p)
        end
    end
end

"""
    κrealisation!(chain::Chain, i::Integer)

Sample the precision parameters for iteration `ì` using Gibbs Sampling.

"""
function κrealisation!(chain::Chain, i::Integer)
    for κp in 1:size(chain.κ, 2)
        p = length(1:κp) - sum(chain.θisRegional[1:κp]) + κp
        α = 1 + (prod(chain.g.G.gridSize) - chain.g.p) / 2
        β = 1 / 10000 +  chain.θ[i, p, :]' * chain.g.G.W * chain.θ[i, p, :] / 2
        chain.κ[i, κp] = rand(Gamma(α, 1 / β ))
    end
end

"""
    fitbayes!(chain::Chain)

Fit the generalized extreme value parameters for every grid cell under the Bayesian paradigm.

"""
function run!(chain::Chain)
    setinitialvalues!(chain)

    @showprogress for i in 2:chain.k
        copylast!(chain, i)

        MHrealisation!(chain, i)

        κrealisation!(chain, i)
    end
end
