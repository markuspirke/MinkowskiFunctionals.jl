using MinkowskiFunctionals
using Distributions
using HDF5
using StatsBase
using DataStructures
using Test

const SAMPLES_DIR = joinpath(@__DIR__, "samples")

@testset "distribution" begin


    # TEST AREA DISTRIBUTION
    n, λ, ρ = 2, 10.0, 10
    dA = AreaDistribution(n^2, λ, ρ)
    @test dA.n == 4
    @test dA.p.p ≈ 0.5420702855281477
    @test dA.p.n == 4
    counts_map = CountsMap(10*ones(2, 2))
    @test abs(compatibility(dA, counts_map)) ≈ 0.1303159919284611
    @test abs(compatibility(dA, counts_map.pixels)) ≈ 0.1303159919284611

    @test 1.0 ≈ sum(pdf(dA))
    @test mean(dA) ≈ dA.n * dA.p.p

    λs = 10.0*ones(2,2)
    dA_nonhomogenous = AreaDistribution(λs, ρ)
    @test 1.0 ≈ sum(pdf(dA_nonhomogenous))
    @test mean(dA) ≈ mean(dA_nonhomogenous)
    @test abs(compatibility(dA_nonhomogenous, counts_map)) ≈ 0.1303159919284611
    @test abs(compatibility(dA_nonhomogenous, counts_map.pixels)) ≈ 0.1303159919284611


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


    # TEST MARGINAL DISTRIBUTIONS
    n = 3
    λ = 10.0
    ρ = 10
    Ω = DensityOfStates(n)
    d = MinkowskiDistribution(Ω, λ, ρ)
    d_A_marginalized = marginalize(d, [:P, :χ])
    d_A = AreaDistribution(n^2, λ, ρ)
    @test pdf(d_A, 0) ≈ pdf(d_A_marginalized, (A=0,))
    @test pdf(d_A, 1) ≈ pdf(d_A_marginalized, (A=1,))
    @test pdf(d_A, 2) ≈ pdf(d_A_marginalized, (A=2,))
    @test pdf(d_A, 3) ≈ pdf(d_A_marginalized, (A=3,))
    @test pdf(d_A, 4) ≈ pdf(d_A_marginalized, (A=4,))
    @test pdf(d_A, 5) ≈ pdf(d_A_marginalized, (A=5,))
    @test pdf(d_A, 6) ≈ pdf(d_A_marginalized, (A=6,))
    @test pdf(d_A, 7) ≈ pdf(d_A_marginalized, (A=7,))
    @test pdf(d_A, 8) ≈ pdf(d_A_marginalized, (A=8,))
    @test pdf(d_A, 9) ≈ pdf(d_A_marginalized, (A=9,))
end
