@enum VarType time rdm

function genVar(t::VarType, n::Integer)::Vector{Float64}
    return t == time ? collect(1:n) : rand(Uniform(), n)
end

function genParam(m::Tuple{<:Integer, <:Integer}, p::Integer, κ::Real, n::Integer,
        ev::Vector{Tuple{VarType, T}} where T<:Real, δ::Real = 0)::Tuple{Array{Real, 2}, Vector{Real}, Array{Real, 2}, Array{Real, 2}}
    F = GMRF.iGMRF(m..., p, κ)
    β₀ = δ .+ GMRF.rand(F)

    βᵢs = zeros(prod(m), 0)
    x = zeros(n, 0)
    for (t, κi) in ev
        βᵢs = hcat(βᵢs, GMRF.rand(GMRF.iGMRF(m..., p, κi)))
        x = hcat(x, genVar(t, n))
    end

    θ = ones(prod(m), n) .* β₀ + βᵢs * x'

    return x, β₀, βᵢs, θ
end

function genGrid(m::Tuple{Integer, Integer}, p::Integer, κμ::Real, κϕ::Real, κξ::Real, n::Integer;
        evμ::Vector{Tuple{VarType, T}} where T<:Real = Vector{Tuple{VarType, Real}}(),
        evϕ::Vector{Tuple{VarType, T}} where T<:Real = Vector{Tuple{VarType, Real}}(),
        evξ::Vector{Tuple{VarType, T}} where T<:Real = Vector{Tuple{VarType, Real}}())::Tuple{BlockMaximaGrid, Array{Real, 2}, Vector{Array{Real, 2}}}

    μx, μβ₀, μβᵢs, μ = genParam(m, p, κμ, n, evμ, 10)

    ϕx, ϕβ₀, ϕβᵢs, ϕ = genParam(m, p, κϕ, n, evϕ)

    ξx, ξβ₀, ξβᵢs, ξ = genParam(m, p, κξ, n, evξ, 0.1)

    gevs = GeneralizedExtremeValue.(μ, exp.(ϕ), ξ)
    y = rand.(gevs)

    return BlockMaximaGrid(m, p, y), hcat(μβ₀, μβᵢs, ϕβ₀, ϕβᵢs, ξβ₀, ξβᵢs), [μx, ϕx, ξx]
end
