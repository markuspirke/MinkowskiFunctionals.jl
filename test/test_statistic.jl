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
    eccdf = ECCDF(λ, L, ecdf(xs), N)

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
    eccdf_loaded = read_eccdf("eccdf_lambda=10.0.h5")
    @test λ == eccdf_loaded.λ
    @test L == eccdf_loaded.L
    @test N == eccdf_loaded.N
    @test xs == eccdf_loaded.ts
    @test eccdf.pvalues == eccdf_loaded.pvalues

    L = 5
    λ = 10.0
    d = Dict(1 => AreaDistribution(L^2, λ, 1))
    dd = DefaultDict(0, d)
    for ρ in 2:100
        dd[ρ] = AreaDistribution(L^2, λ, ρ)
    end
    ts = [calc_ts(dd, CountsMap(L, λ)) for _ in 1:100]
    e_cdf = ecdf(ts)
    eccdf = ECCDF(λ, L, e_cdf, 100)
    # @test 1 - e_cdf(1.0) ≈ eccdf(1.0) these will not work, as e_cdf has no equidistance bins, think about what to do here...
    # @test 1 - e_cdf(20.0) ≈ eccdf(20.0)
    counts_map = CountsMap(16, λ)
    update_distributions!(dd, counts_map)
    mink_map = MinkowskiMap(counts_map, dd, eccdf)
    @test !isnan(maximum(p2σ(mink_map)))
    mink_map = MinkowskiMap(counts_map, dd, e_cdf)
    @test !isnan(maximum(p2σ(mink_map)))

    counts_map = CountsMap(5, 30.0)
    mink_map = MinkowskiMap(counts_map, dd, e_cdf)
    @test 0.01 ≈ mink_map[1, 1]
    mink_map = MinkowskiMap(counts_map, dd, eccdf)
    @test 0.01 ≈ mink_map[1, 1]


end
