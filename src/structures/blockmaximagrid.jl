struct BlockMaximaGrid
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

"""
    Base.show(io::IO, obj::BlockMaximaGrid)

Displays an object of type BlockMaximaGrid.

"""
function Base.show(io::IO, obj::BlockMaximaGrid)
    showGrid(io, obj)
end

"""
    showGrid(io::IO, obj::BlockMaximaGrid; prefix::String = "")

Displays an object of type BlockMaximaGrid with a prefix.

"""
function showGrid(io::IO, obj::BlockMaximaGrid; prefix::String = "")
    println(io, prefix, "BlockMaximaGrid")
    println(io, prefix, "p :\t", obj.p)
    println(io, prefix, "G :")
    GMRF.showGridStructure(io, obj.G, prefix = prefix * "\t" )
    println(io)
    println(io, prefix, "y :\t", typeof(obj.y), "[", size(obj.y, 1), ", ", size(obj.y, 2), "]")
end
