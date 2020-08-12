@testset "util.jl" begin
    @testset "multiply(m, v)" begin
        # 1x1 matrix times 1x1 vector
        @test RegionalExtremes.multiply(ones(1, 1), ones(1)) == [1]

        # 2x2 matrix times 2x1 vector
        @test RegionalExtremes.multiply([1 2; 3 4], [5, 6]) == [17, 39]

        # 1x3 matrix times 3x1 vector
        @test RegionalExtremes.multiply([7 8 9], [10, 11, 12]) == [266]
    end
end
