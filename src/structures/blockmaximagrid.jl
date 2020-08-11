struct BlockMaximaGrid # TODO : println
    p::Integer
    G::GMRF.GridStructure
    y::Array{<:Real, 2}
end

Base.Broadcast.broadcastable(obj::BlockMaximaGrid) = Ref(obj)

"""
    BlockMaximaGrid(m::Tuple{Integer, Integer}, p::Integer, y::Array{<:Real, 2})::BlockMaximaGrid

Build a BlockMaximaGrid.

"""
function BlockMaximaGrid(m::Tuple{Integer, Integer}, p::Integer, y::Array{<:Real, 2})::BlockMaximaGrid
    @assert prod(m) == size(y, 1) "Number of grid cells must correspond to the number of rows in the data matrix."

    return BlockMaximaGrid(p, iGMRF(m..., p, 1).G, y)
end
