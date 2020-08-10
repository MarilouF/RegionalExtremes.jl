function setinitialvalues!(chain::Chain)
    gumbels = Extremes.gumbelfitpwm.([chain.g.y[i, :] for i in 1:size(chain.g.y, 1)])
    for i in 1:length(gumbels)
        for p in 1:3
            chain.θ[1, chain.index[p][1], i] = gumbels[i].θ̂[p]
        end
    end
    for p in 1:size(chain.θ, 2)
        chain.θlogpdf[1, p, :] = local_loglikelihood(chain, chain.θ[1, :, :], collect(1:prod(chain.g.G.gridSize)), 1)
        if chain.θisRegional[p]
            κ = chain.θ[1, p, :]' * chain.g.G.W * chain.θ[1, p, :] ./ (prod(chain.g.G.gridSize) - 1)
            chain.κ[1, sum(chain.θisRegional[1:p])] = κ == 0 ? 1.0 : κ # TODO : 1.0 ?

            chain.θlogpdf[1, p, :] .+= regional_loglikelihood(chain, chain.θ[1, :, :], collect(1:prod(chain.g.G.gridSize)), 1, p)
        end
    end
end

function copylast!(chain::Chain, i::Integer)
    chain.θ[i, :, :] = chain.θ[i - 1, :, :]
    chain.κ[i, :] = chain.κ[i - 1, :]
    chain.θlogpdf[i, :, :] = chain.θlogpdf[i - 1, :, :]
end

function candidates(chain::Chain, i::Integer, subset::Vector{<:Integer}, param::Integer)::Array{<:Real, 2}
    c = copy(chain.θ[i, :, :])
    c[param, subset] = rand.(Normal.(chain.θ[i, param, subset], chain.θstep[param]))
    return c
end

function MHupdatechain!(chain::Chain, i::Integer, subset::Vector{<:Integer}, param::Integer, c::Array{<:Real, 2}, logpdfP::Vector{<:Real})
    logα = min.(0, logpdfP - chain.θlogpdf[i, param, subset])
    logu = log(rand(Uniform()))

    a = logu .< logα

    chain.θacc[i, param, subset] = a
    chain.θ[i, param, subset[a]] = c[param, subset[a]]
    chain.θlogpdf[i, param, subset[a]] = logpdfP[a]
end

function MHlocal!(chain::Chain, i::Integer, p::Integer)
    set = collect(1:size(chain.g.y, 1))
    c = candidates(chain, i, set, p)
    logpdf = local_loglikelihood(chain, c, set, i)
    MHupdatechain!(chain, i, set, p, c, logpdf)
end

function MHregional!(chain::Chain, i::Integer, p::Integer)
    for subset in chain.g.G.condIndSubset
        c = candidates(chain, i, subset, p)
        logpdf = local_loglikelihood(chain, c, subset, i) .+ regional_loglikelihood(chain, c, subset, i, p)
        MHupdatechain!(chain, i, subset, p, c, logpdf)
    end
end

function MHrealisation!(chain::Chain, i::Integer)
    for p in 1:size(chain.θ, 2)
        if chain.θisRegional[p]
            MHregional!(chain, i, p)
        else
            MHlocal!(chain, i, p)
        end
    end
end

function κrealisation!(chain, i::Integer)
    for κp in 1:size(chain.κ, 2)
        p = length(1:κp) - sum(chain.θisRegional[1:κp]) + κp
        α = 1 + (prod(chain.g.G.gridSize) - chain.g.p) / 2
        β = 1 / 10000 +  chain.θ[i, p, :]' * chain.g.G.W * chain.θ[i, p, :] / 2
        chain.κ[i, κp] = rand(Gamma(α, 1 / β ))
    end
end

function fitbayes!(chain::Chain)
    setinitialvalues!(chain)

    for i in 2:chain.k
        copylast!(chain, i)

        MHrealisation!(chain, i)

        κrealisation!(chain, i)
    end
end
