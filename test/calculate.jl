using Test
using MinkowskiFunctionals

@testset "calculate" begin
    image = [1 1 1; 1 0 1; 1 1 1]
    M = MinkowskiFunctional(BWMap(image))

    @test M.A == 8
    @test M.P == 16
    @test M.χ == 0

    M = MinkowskiFunctional(image)

    @test M.A == 8
    @test M.P == 16
    @test M.χ == 0

    image = zeros(42, 42)
    image[42] = 1
    M = MinkowskiFunctional(BWMap(image))

    @test M.A == 1
    @test M.P == 4
    @test M.χ == 1

end
