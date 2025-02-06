using Test
using MinkowskiFunctionals
using Distributions

const SAMPLES_DIR = joinpath(@__DIR__, "samples")

@testset "minkowski_map" begin
    # ρs = 9:10
    # λ = 10
    # n = 3
    # distribution = [PoissonMinkowskiDistributions(n, λ, ρ) for ρ in ρs]
    # cmap = CountsMap([10 10 10; 10 0 10; 10 10 10])
    # mmap = MinkowskiMap(cmap, distribution, :A)
    # @test 3 == mmap.L
    # @test 10 == mmap.λ
    # @test 9:10 == mmap.ρ
    # @test (1, 1) == size(mmap.pixels)
    # @test 1.3459674701327748 ≈ mmap.pixels[1]

    # mmap = MinkowskiMap(cmap, distribution, :P)
    # @test 0.4658702309200946 ≈ mmap.pixels[1]

    # cmap = CountsMap(
    #     [10 10 10 10 10 ;
    #      10 10 10 10 10 ;
    #      10 10 0 10 10 ;
    #      10 10 10 10 10 ;
    #      10 10 10 10 10
    #     ]
    #    )
    # mmap = MinkowskiMap(cmap, distribution, :A)
    # @test  (mmap.pixels .≈ 1.3459674701327748) == ones(3, 3)
    # mmap = MinkowskiMap(cmap, distribution, :P)
    # correct = ones(3, 3) * -9.64327466553287e-17
    # correct[1, 2] = 0.19514784000984248
    # correct[2, 1] = 0.19514784000984248
    # correct[2, 3] = 0.19514784000984248
    # correct[3, 2] = 0.19514784000984248
    # correct[2, 2] = 0.4658702309200946
    # @test correct ≈ mmap.pixels

    Ω = DensityOfStates(joinpath(SAMPLES_DIR, "structure_5x5"))
    h0s = [MinkowskiDistribution(Ω, 10, ρ) for ρ in 10:11]
    counts_map = rand(CountsMap, Poisson(10), 8)
    mink_map = MinkowskiMap(counts_map, h0s)

    mink_map = MinkowskiMap(1, 1, 1, 0.5 * σ2p(5) .* ones(2, 2))
    mink_map = correct_trials(mink_map, 2)
    @test mink_map.pixels ≈ σ2p(5) .* ones(2, 2)
    mink_map = p2σ(mink_map)
    @test mink_map.pixels ≈ 5.0 .* ones(2, 2)
end
