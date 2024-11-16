using Test
using MinkowskiFunctionals

@testset "density_of_states" begin
    Ω = DensityOfStates(2)
    @test 1 == Ω[MinkowskiFunctional(0, 0, 0)]
    @test 4 == Ω[MinkowskiFunctional(1, 4, 1)]
    @test 4 == Ω[MinkowskiFunctional(2, 6, 1)]
    @test 2 == Ω[MinkowskiFunctional(2, 8, 1)]
    @test 4 == Ω[MinkowskiFunctional(3, 8, 1)]
    @test 1 == Ω[MinkowskiFunctional(4, 8, 1)]
    # @test 1 == Ω[MinkowskiFunctional(0, 0, 0)]


end
