using Test
using MinkowskiFunctionals

@testset "minkowski_map" begin
    @test 1.0 ≈ λ_lima(1, 1)
    @test 1.0 ≈ λ_lima(10, 10)
    @test (1/2 * 3)^1 * (1/2 * 3/2)^2 ≈ λ_lima(1, 2)
    @test 0.0 ≈ significance_lima(1, 1)
    @test 0.0 ≈ significance_lima(10, 10)
    @test sqrt(-2log((1/2 * 3)^1 * (1/2 * 3/2)^2)) ≈ significance_lima(1, 2)

    counts_map = CountsMap(Int.(ones(9, 9)))
    λ = 1
    L = 3
    @test zeros(7, 7) ≈ lima_map(counts_map, λ, L)
    @test zeros(7, 7) ≈ lima_map_roundkernel(counts_map, λ, L)
end
