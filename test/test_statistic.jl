using MinkowskiFunctionals
using Distributions
using StatsBase
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
    eccdf = ECCDF(λ, L, ecdf(xs), xs)

    @test 0.9 ≈ eccdf(1.0)
    @test 0.8 ≈ eccdf(2.0)
    @test 0.7 ≈ eccdf(3.0)
    @test 0.6 ≈ eccdf(4.0)
    @test 0.5 ≈ eccdf(5.0)
    @test 0.4 ≈ eccdf(6.0)
    @test 0.3 ≈ eccdf(7.0)
    @test 0.2 ≈ eccdf(8.0)
    @test 0.1 ≈ eccdf(9.0)
    @test 0.0 ≈ eccdf(10.0)
end
