@testset "blockmaximagrid.jl" begin
    @testset "BlockMaximaGrid(m, p, y)" begin
        # Grid size not corresponding to the number of data rows throws
        @test_throws AssertionError BlockMaximaGrid((1, 2), 1, [1 0])

        # Grid of order 1
        m = (1, 2)
        p = 1
        y = [1 0; 1 0]
        grid = BlockMaximaGrid(m, p, y)

        @test grid.G.gridSize == m
        @test grid.p == p
        @test grid.y == y

        # Grid of order 2
        m = (1, 2)
        p = 2
        y = [1 0; 1 0]
        grid = BlockMaximaGrid(m, p, y)

        @test grid.G.gridSize == m
        @test grid.p == p
        @test grid.y == y
    end
    g = BlockMaximaGrid((1, 2), 1, [1 2 3; 4 5 6])
    @testset "Base.show(io, obj)" begin
        # Print does not throw
        io = IOBuffer()
        @test_logs Base.show(io, g)
    end
    @testset "showGrid(io, obj; prefix)" begin
        # Print does not throw
        io = IOBuffer()
        @test_logs RegionalExtremes.showGrid(io, g, prefix = "\t")
    end
end
