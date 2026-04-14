using MinkowskiFunctionals
using Distributions
using StatsBase
using DataStructures
using HDF5
using Test

const SAMPLES_DIR = joinpath(@__DIR__, "samples")

@testset "test-statistic" begin
    λ = 10.0
    L = 15
    xs = 1:10
    N = 10
    eccdf = ECCDF(λ, L, N, ecdf(collect(xs)))

    @test 0.9 ≈ eccdf(1.0)
    @test 0.8 ≈ eccdf(2.0)
    @test 0.7 ≈ eccdf(3.0)
    @test 0.6 ≈ eccdf(4.0)
    @test 0.5 ≈ eccdf(5.0)
    @test 0.4 ≈ eccdf(6.0)
    @test 0.3 ≈ eccdf(7.0)
    @test 0.2 ≈ eccdf(8.0)
    @test 0.1 ≈ eccdf(9.0)
    @test 0.1 ≈ eccdf(10.0)
    @test 0.1 ≈ eccdf(11.0)

    write_eccdf(".", eccdf)
    eccdf_loaded = read_eccdf("eccdf_lambda=$(λ)_L=$(L).h5")
    @test λ == eccdf_loaded.λ
    @test L == eccdf_loaded.L
    @test N == eccdf_loaded.N
    @test 0.6 ≈ eccdf_loaded(4.0)
    @test 0.5 ≈ eccdf_loaded(5.0)

    λ = 10.0
    counts_map = CountsMap(16, λ)
    Ω = DensityOfStates(3)
    lut = MinkowskiPValueLookup(Ω)
    eccdf = ECCDF(λ, lut, 100)

    mink_map = MinkowskiMap(counts_map, λ, lut, eccdf)
    @test !isnan(maximum(p2σ(mink_map)))
    background = Background(λ * ones(16, 16))
    mink_map = MinkowskiMap(counts_map, background, lut, 100)
    @test !isnan(maximum(p2σ(mink_map)))
end
