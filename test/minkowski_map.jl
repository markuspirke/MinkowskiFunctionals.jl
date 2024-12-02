
using Test
using MinkowskiFunctionals

@testset "minkowski_map" begin
    ρs = 9:10
    λ = 10
    n = 3
    distribution = [PoissonMinkowskiDistributions(n, λ, ρ) for ρ in ρs]
    cmap = CountsMap([10 10 10; 10 0 10; 10 10 10])
    mmap = MinkowskiMap(cmap, distribution, :A)
    @test 3 == mmap.L
    @test 10 == mmap.λ
    @test 9:10 == mmap.ρ
    @test (1, 1) == size(mmap.pixels)
    @test 1.3459674701327748 ≈ mmap.pixels[1]

    mmap = MinkowskiMap(cmap, distribution, :P)
    @test 0.4658702309200946 ≈ mmap.pixels[1]

    cmap = CountsMap(
        [10 10 10 10 10 ;
         10 10 10 10 10 ;
         10 10 0 10 10 ;
         10 10 10 10 10 ;
         10 10 10 10 10
        ]
       )
    mmap = MinkowskiMap(cmap, distribution, :A)
    @test  (mmap.pixels .≈ 1.3459674701327748) == ones(3, 3)
    mmap = MinkowskiMap(cmap, distribution, :P)
    correct = ones(3, 3) * -9.64327466553287e-17
    correct[1, 2] = 0.19514784000984248
    correct[2, 1] = 0.19514784000984248
    correct[2, 3] = 0.19514784000984248
    correct[3, 2] = 0.19514784000984248
    correct[2, 2] = 0.4658702309200946
    @test correct ≈ mmap.pixels
end
