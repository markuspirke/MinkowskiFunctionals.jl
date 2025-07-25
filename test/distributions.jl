using MinkowskiFunctionals
using Distributions
using HDF5
using StatsBase
using DataStructures
using Test

const SAMPLES_DIR = joinpath(@__DIR__, "samples")

@testset "distribution" begin


    # TEST AREA DISTRIBUTION
    n, λ, ρ = 2, 10, 10
    dA = AreaDistribution(n^2, λ, ρ)
    @test dA.n == 4
    @test dA.p.p ≈ 0.5420702855281477
    @test dA.p.n == 4
    counts_map = CountsMap(10*ones(2, 2))
    @test abs(compatibility(dA, counts_map)) ≈ 0.1303159919284611
    @test abs(compatibility(dA, counts_map.pixels)) ≈ 0.1303159919284611

    @test 1.0 ≈ sum(pdf(dA))


    n, λ, ρ = 3, 10, 10
    Ω = DensityOfStates(n)
    @test 2^(n^2) == sum(values(Ω.data))
    P = MinkowskiDistribution(Ω, λ, ρ)
    @test 1.0 ≈ sum(pdf(P))



    # TEST IO
    n = 3
    λ = 10.0
    ρ = 10
    Ω = DensityOfStates(n)
    d = MinkowskiDistribution(Ω, λ, ρ, pvalues=false)
    h5open("foo.h5", "w") do h5f
        append!(h5f, d)
    end

    h5open("foo.h5", "r") do h5f
        d_loaded = MinkowskiDistribution(h5f, n, λ, ρ)
        @test d_loaded.n == n
        @test d_loaded.ρ == ρ
        @test d_loaded.λ == λ
        @test d_loaded.p == d.p
        @test ismissing(d_loaded.pvalues)
    end

    d = MinkowskiDistribution(Ω, λ, ρ, pvalues=true)
    h5open("foo.h5", "w") do h5f
        append!(h5f, d)
    end

    h5open("foo.h5", "r") do h5f
        d_loaded = MinkowskiDistribution(h5f, n, λ, ρ)
        @test d_loaded.n == n
        @test d_loaded.ρ == ρ
        @test d_loaded.λ == λ
        @test d_loaded.p == d.p
        d_pvalues = sort(collect(values(d.pvalues)))
        d_loaded_pvalues = sort(collect(values(d_loaded.pvalues)))
        @test d_loaded_pvalues ≈ d_pvalues
    end

    rm("foo.h5")

    Ω = DensityOfStates(joinpath(SAMPLES_DIR, "structure_5x5"))
    d = MinkowskiDistribution(Ω, 10.0, 10)
    write_pvalues(".", d)
    pvalues = read_pvalues("lambda=10.0_rho=10.dat")
    @test pvalues == d.pvalues
    rm("lambda=10.0_rho=10.dat")
    d_pvalues = Dict(d.pvalues)

    x = [
    15  12   8  12  20;
    16   8   8  12   6;
    12   8  16  10  12;
    12  10  13  14   7;
    16   9  16   8  13
    ]
    x_counts_map = CountsMap(x)

    @test 0.15485799310649023 ≈ compatibility(d_pvalues, 10, x)
    @test 0.15485799310649023 ≈ compatibility(d_pvalues, 10, x_counts_map)



    # TEST SIGN CONVENTION
    # n, λ, ρ = 9, 10, 11
    # counts_map = CountsMap([12 12 12; 12 9 9; 12 12 12])
    # dA = AreaDistribution(n, λ, ρ)
    # @test MinkowskiFunctionals.get_sign(counts_map[:, :], dA) == 1.0
    # counts_map = CountsMap([12 9 9; 12 9 9; 9 9 12])
    # @test MinkowskiFunctionals.get_sign(counts_map[:, :], dA) == -1.0



    # MINKOWSKI DEVIATION STRENGTH
    # p_counter = Accumulator{MinkowskiFunctional, Float64}()
    # p_counter[MinkowskiFunctional(0,0,0)] = 0.15865525393145696
    # p_counter[MinkowskiFunctional(1,0,0)] = 0.15865525393145696
    # p_counter[MinkowskiFunctional(2,0,0)] = 0.6826894921370861

    # D = MinkowskiDistribution(1, 1, 1, p_counter)
    # @test p2σ(compatibility(D, MinkowskiFunctional(0,0,0))) ≈ 1.0 # in sigma
    # @test p2σ(compatibility(D, MinkowskiFunctional(1,0,0))) ≈ 1.0 # in sigma
    # @test p2σ(compatibility(D, MinkowskiFunctional(2,0,0))) ≈ 0.0 # in sigma

    # h5open(joinpath(SAMPLES_DIR, "test_lookuptable.h5"), "w") do h5f
    #     append!(h5f, D)
    # end

    # D_loaded = MinkowskiDistribution(joinpath(SAMPLES_DIR, "test_lookuptable.h5"), 1, 1, 1)
    # @test collect(keys(D.p)) == collect(keys(D_loaded.p))
    # @test collect(keys(D.α)) == collect(keys(D_loaded.α))
    # @test D.n == D_loaded.n
    # @test D.λ == D_loaded.λ
    # @test D.ρ == D_loaded.ρ
    # @test D.p == D_loaded.p
    # @test D.α == D_loaded.α

    # # CHECK WHETHER APPENDING ALSO WORKS
    # D = MinkowskiDistribution(2, 1, 1, p_counter)

    # h5open(joinpath(SAMPLES_DIR, "test_lookuptable.h5"), "cw") do h5f
    #     append!(h5f, D)
    #     @test parse.(Int, keys(h5f)) == [1, 2]
    # end

    # D_loaded = MinkowskiDistribution(joinpath(SAMPLES_DIR, "test_lookuptable.h5"), 2, 1, 1)
    # @test collect(keys(D.p)) == collect(keys(D_loaded.p))
    # @test collect(keys(D.α)) == collect(keys(D_loaded.α))
    # @test D.n == D_loaded.n
    # @test D.λ == D_loaded.λ
    # @test D.ρ == D_loaded.ρ
    # @test D.p == D_loaded.p
    # @test D.α == D_loaded.α

    # D = MinkowskiDistribution(2, 2, 1, p_counter)
    # h5open(joinpath(SAMPLES_DIR, "test_lookuptable.h5"), "cw") do h5f
    #     append!(h5f, D)
    #     @test parse.(Int, keys(h5f)) == [1, 2]
    # end
    # D_loaded = MinkowskiDistribution(joinpath(SAMPLES_DIR, "test_lookuptable.h5"), 2, 2, 1)
    # @test collect(keys(D.p)) == collect(keys(D_loaded.p))
    # @test collect(keys(D.α)) == collect(keys(D_loaded.α))
    # @test D.n == D_loaded.n
    # @test D.λ == D_loaded.λ
    # @test D.ρ == D_loaded.ρ
    # @test D.p == D_loaded.p
    # @test D.α == D_loaded.α

    # D = MinkowskiDistribution(2, 2, 2, p_counter)
    # h5open(joinpath(SAMPLES_DIR, "test_lookuptable.h5"), "cw") do h5f
    #     append!(h5f, D)
    #     @test parse.(Int, keys(h5f)) == [1, 2]
    # end
    # D_loaded = MinkowskiDistribution(joinpath(SAMPLES_DIR, "test_lookuptable.h5"), 2, 2, 2)
    # @test collect(keys(D.p)) == collect(keys(D_loaded.p))
    # @test collect(keys(D.α)) == collect(keys(D_loaded.α))
    # @test D.n == D_loaded.n
    # @test D.λ == D_loaded.λ
    # @test D.ρ == D_loaded.ρ
    # @test D.p == D_loaded.p
    # @test D.α == D_loaded.α

    #TEST MARGINAL DISTRIBUTIONS
    # p_counter = Accumulator{MinkowskiFunctional, Float64}()
    # p_counter[MinkowskiFunctional(0,0,0)] = 0.1
    # p_counter[MinkowskiFunctional(0,1,0)] = 0.2
    # p_counter[MinkowskiFunctional(1,0,0)] = 0.3
    # p_counter[MinkowskiFunctional(1,1,0)] = 0.4
    # D = MinkowskiDistribution(1, 1, 1, p_counter)

    # @test 0.1 ≈ D.α[MinkowskiFunctional(0,0,0)]
    # @test 0.3 ≈ D.α[MinkowskiFunctional(0,1,0)]
    # @test 0.6 ≈ D.α[MinkowskiFunctional(1,0,0)]
    # @test 1.0 ≈ D.α[MinkowskiFunctional(1,1,0)]

    # D_A = marginalize(D, (:P, :χ))
    # @test 0.3 ≈ D_A.p[(A =0,)]
    # @test 0.7 ≈ D_A.p[(A =1,)]
    # @test 0.3 ≈ D_A.α[(A =0,)]
    # @test 1.0 ≈ D_A.α[(A =1,)]

    # Ω = DensityOfStates(joinpath(SAMPLES_DIR, "structure_5x5"))
    # D = MinkowskiDistribution(Ω, 10, 10)
    # d_poisson = Distributions.Poisson(D.λ)
    # p = 1 -cdf(d_poisson, D.ρ-1)
    # D_A_binomial = Binomial(D.n^2, p)
    # D_A = marginalize(D, (:P, :χ))

    # ps_equal = Bool[]
    # for x in 1:25
    #     push!(ps_equal, pdf(D_A_binomial, x) ≈ D_A.p[(A=x,)])
    # end
    # @test sum(ps_equal) == 25

    # ps = pdf(D_A_binomial)
    # αs_equal = Bool[]
    # for x in 1:25
    #     push!(αs_equal, round(sum(ps[ps .<= pdf(D_A_binomial, x)]), digits=10) ≈ round(D_A.α[(A=x,)], digits=10))
    # end
    # @test sum(αs_equal) == 25
end
