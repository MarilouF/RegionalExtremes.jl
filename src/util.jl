# Simple matrix and vector multiplication, faster than using *
function multiply(m::Array{<:Real, 2}, v::Vector{<:Real})::Vector{<:Real}
    sol = m[:, 1] * v[1]
    for i in 2:size(m, 2)
        sol += m[:, i] * v[i]
    end
    return sol
end
