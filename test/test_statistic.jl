using MinkowskiFunctionals
using Distributions
using StatsBase
using DataStructures
using Test

const SAMPLES_DIR = joinpath(@__DIR__, "samples")

@testset "test-statistic" begin
    x = AreaDistributionX(25, 10.0, 5)
    y = AreaDistribution(25, 10.0, 5)
    countsmap = CountsMap(5, 10.0)
    @test compatibility(x, countsmap) ≈ compatibility(y, countsmap)

    x = AreaDistributionX(25, 10.0, 10)
    y = AreaDistribution(25, 10.0, 10)
    countsmap = CountsMap(5, 10.0)
    @test compatibility(x, countsmap) ≈ compatibility(y, countsmap)

    x = AreaDistributionX(25, 10.0, 15)
    y = AreaDistribution(25, 10.0, 15)
    countsmap = CountsMap(5, 10.0)
    @test compatibility(x, countsmap) ≈ compatibility(y, countsmap)

    f = h5open("foo.h5", "w")
    append!(f, x)
    y = AreaDistributionX(f, x.ρ)
    @test x.n == y.n
    @test x.λ == y.λ
    @test x.p == y.p
    @test x.pvalues == y.pvalues
    close(f)
    rm("foo.h5")

    λ = 10.0
    L = 15
    xs = 1:10
    eccdf = ECCDF(λ, L, ecdf(xs), 10)

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

    L = 5
    λ = 10.0
    d = Dict(1 => AreaDistributionX(L^2, λ, 1))
    dd = DefaultDict(0, d)
    ts = [calc_ts!(dd, CountsMap(L, λ)) for _ in 1:100]
    e_cdf = ecdf(ts)
    eccdf = ECCDF(λ, L, e_cdf, 100)
    # @test 1 - e_cdf(1.0) ≈ eccdf(1.0) these will not work, as e_cdf has no equidistance bins, think about what to do here...
    # @test 1 - e_cdf(20.0) ≈ eccdf(20.0)
    counts_map = CountsMap(16, λ)
    update_distributions!(dd, counts_map)
    mink_map = MinkowskiMap(counts_map, dd, eccdf)
    @test !isnan(maximum(p2σ(mink_map)))

end
