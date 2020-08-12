# RegionalExtremes.jl

## Dependency
The module GMRF.jl is an unpublished dependency of this project.
It must be downloaded alongside this module in order to use the latter.

## Tutorial
To start this tutorial, import the following dependencies:

    using RegionalExtremes # Current library
    using Plots, Statistics # Displaying of results

### Stationary Model
The first (and simplest) option will be to generate and estimate a stationary model.\
The alternative would be to use a non-stationary model which will be presented in the next section.

We'll generate a 10 x 15 grid with 300 years of maxima recorded.\
This grid of first order will have a precision parameter for each generalized extreme value distribution parameter.\
These parameters will be 50, 50 and 3000 for the location, scale and shape parameters respectively.

    gridSize = (10, 15)
    ndata = 300
    order = 1
    κs = [50.0, 50.0, 3000.0]

The grid will be generated using the `genGrid` function in the `RegionalExtremes` library.\
The `g` is the grid with its generated maxima.\
The `θ` is the value of the generated generalized extreme value parameters for every grid cell.\
The `x` is the generated explanatory variables (none here).

    g, θ, x = genGrid(gridSize, order, κs..., ndata)

The Gibbs sampling will be run for 3000 iterations (500 of burn in).\
The step for the random walk for the generalized extreme value parameter (location, scale and shape) will be 0.2, 0.15 and 0.04.\
By default, all the estimated parameters are considered in the regional model. It is possible to indicate that a parameter is to be estimated using the local data only. Here, we will do so for the scale parameter.

    niteration = 3000
    burnin = 500
    steps = [0.2, 0.15, 0.04]
    regional = [true, false, true]

In order to estimate the parameters, the `fit` function is called.

    c = fit(gridSize, order, g.y, niteration, burnin, steps, regional = regional)

The value obtained, `c`, is an MCMC chain resulting from the Gibbs sampling algorithm.\
Visualisation of the results yielded from this chain will be presented in the **Using the Resulting Chain** section.

### Non-stationary Model
The second option is to generate and estimate a non-stationary model.


We'll generate a 10 x 15 grid with 300 years of non-stationary maxima recorded.\
This grid will be of first order.

The location parameter (of the generalized extreme value distribution) of each grid cell will vary according to a certain explanatory variable.\
This explanatory variable will be randomly sampled from a uniform distribution [0, 1].\
Therefore, two parameters will be generated for the location parameter of every grid cell.\
The first one will have a precision parameter of 50 and the second one (associated to the explanatory variable) will have a precision parameter of 10.

The scale and shape parameters will be stationary. They will have a precision parameter of 50 and 3000 respectively.

    gridSize = (10, 15)
    ndata = 300
    order = 1
    κs = [50.0, 10.0, 50.0, 3000.0]
    explVar = (RegionalExtremes.rdm, κs[2]) # rdm indicates a uniform distribution for the generator

The grid will be generated using the `genGrid` function in the `RegionalExtremes` library.\
The `g` is the grid with its generated maxima.\
The `θ` is the value of the generated generalized extreme value parameters for every grid cell.\
The `x` is the generated explanatory variables (The uniformly distributed variable).

    g, θ, x = genGrid(gridSize, order, κs[1], κs[3], κs[4], ndata, evμ = [explVar])

The Gibbs sampling will be run for 3000 iterations (500 of burn in).\
The step for the random walk for the generalized extreme value parameter (location1, location2, scale and shape) will be 0.2, 0.3, 0.15 and 0.04.\
By default, all the estimated parameters are considered in the regional model. It is possible to indicate that a parameter is to be estimated using the local data only. Here, we will do so for the shape parameter.

    niteration = 3000
    burnin = 500
    steps = [0.2, 0.3, 0.15, 0.15]
    regional = [true, true, true, false]

In order to estimate the parameters, the `fit` function is called.

    c = fit(gridSize, order, g.y, niteration, burnin, steps, regional = regional, μcov = x[1])

The value obtained, `c`, is an MCMC chain resulting from the Gibbs sampling algorithm.\
Visualisation of the results yielded from this chain will be presented in the next section.

### Using the Resulting Chain
Whether the chain you are using was created from a stationary or non-stationary model, the following code snippets will allow you to visualize your results.

The following code will allow you to visualize the actual vs estimated generalized extreme value parameters with heatmaps.

    np = sum([length(v) for v in c.index])
    plots = []
    for i in 1:np
        push!(plots, heatmap(reshape(θ[:, i], c.g.G.gridSize...)))
        push!(plots, heatmap(reshape(mean(c.θ[c.burnin:end, i, :], dims=1), c.g.G.gridSize...)))
    end
    plot(plots..., layout = (np,2), size = (600, 1000))

The following code will allow you to visualize the actual vs estimated confidence interval for the precision parameters.

    for p in 1:size(c.κ, 2)
        q1 = quantile!(c.κ[c.burnin:end, p], 0.025)
        q2 = quantile!(c.κ[c.burnin:end, p], 0.975)
        println(κs[p - sum(c.θisRegional[1:p]) + p], " - [", q1, ", ", q2, "]")
    end
