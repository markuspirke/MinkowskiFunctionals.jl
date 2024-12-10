using MinkowskiFunctionals
using Distributions
using Test

const SAMPLES_DIR = joinpath(@__DIR__, "samples")

@testset "distribution" begin

    n, λ, ρ = 3, 10, 10
    Ω = DensityOfStates(n)
    @test 2^(n^2) == sum(values(Ω.data))
    P = MinkowskiDistribution(Ω, λ, ρ)
    @test 1.0 ≈ sum(pdf(P))


    # MINKOWSKI DEVIATION STRENGTH
    p_counter = Accumulator{MinkowskiFunctional, Float64}()
    p_counter[MinkowskiFunctional(0,0,0)] = 0.15865525393145696
    p_counter[MinkowskiFunctional(1,0,0)] = 0.15865525393145696
    p_counter[MinkowskiFunctional(2,0,0)] = 0.6826894921370861

    D = MinkowskiDistribution(1, 1, 1, p_counter)
    @test p2σ(compatibility(D, MinkowskiFunctional(0,0,0))) ≈ 1.0 # in sigma
    @test p2σ(compatibility(D, MinkowskiFunctional(1,0,0))) ≈ 1.0 # in sigma
    @test p2σ(compatibility(D, MinkowskiFunctional(2,0,0))) ≈ 0.0 # in sigma

    h5open(joinpath(SAMPLES_DIR, "test_lookuptable.h5"), "w") do h5f
        append!(h5f, D)
    end

    D_loaded = MinkowskiDistribution(joinpath(SAMPLES_DIR, "test_lookuptable.h5"), 1, 1, 1)
    @test collect(keys(D.p)) == collect(keys(D_loaded.p))
    @test collect(keys(D.α)) == collect(keys(D_loaded.α))
    @test D.n == D_loaded.n
    @test D.λ == D_loaded.λ
    @test D.ρ == D_loaded.ρ
    @test D.p == D_loaded.p
    @test D.α == D_loaded.α

    # CHECK WHETHER APPENDING ALSO WORKS
    D = MinkowskiDistribution(2, 1, 1, p_counter)

    h5open(joinpath(SAMPLES_DIR, "test_lookuptable.h5"), "cw") do h5f
        append!(h5f, D)
        @test parse.(Int, keys(h5f)) == [1, 2]
    end

    D_loaded = MinkowskiDistribution(joinpath(SAMPLES_DIR, "test_lookuptable.h5"), 2, 1, 1)
    @test collect(keys(D.p)) == collect(keys(D_loaded.p))
    @test collect(keys(D.α)) == collect(keys(D_loaded.α))
    @test D.n == D_loaded.n
    @test D.λ == D_loaded.λ
    @test D.ρ == D_loaded.ρ
    @test D.p == D_loaded.p
    @test D.α == D_loaded.α

    D = MinkowskiDistribution(2, 2, 1, p_counter)
    h5open(joinpath(SAMPLES_DIR, "test_lookuptable.h5"), "cw") do h5f
        append!(h5f, D)
        @test parse.(Int, keys(h5f)) == [1, 2]
    end
    D_loaded = MinkowskiDistribution(joinpath(SAMPLES_DIR, "test_lookuptable.h5"), 2, 2, 1)
    @test collect(keys(D.p)) == collect(keys(D_loaded.p))
    @test collect(keys(D.α)) == collect(keys(D_loaded.α))
    @test D.n == D_loaded.n
    @test D.λ == D_loaded.λ
    @test D.ρ == D_loaded.ρ
    @test D.p == D_loaded.p
    @test D.α == D_loaded.α

    D = MinkowskiDistribution(2, 2, 2, p_counter)
    h5open(joinpath(SAMPLES_DIR, "test_lookuptable.h5"), "cw") do h5f
        append!(h5f, D)
        @test parse.(Int, keys(h5f)) == [1, 2]
    end
    D_loaded = MinkowskiDistribution(joinpath(SAMPLES_DIR, "test_lookuptable.h5"), 2, 2, 2)
    @test collect(keys(D.p)) == collect(keys(D_loaded.p))
    @test collect(keys(D.α)) == collect(keys(D_loaded.α))
    @test D.n == D_loaded.n
    @test D.λ == D_loaded.λ
    @test D.ρ == D_loaded.ρ
    @test D.p == D_loaded.p
    @test D.α == D_loaded.α

    #TEST MARGINAL DISTRIBUTIONS
    p_counter = Accumulator{MinkowskiFunctional, Float64}()
    p_counter[MinkowskiFunctional(0,0,0)] = 0.1
    p_counter[MinkowskiFunctional(0,1,0)] = 0.2
    p_counter[MinkowskiFunctional(1,0,0)] = 0.3
    p_counter[MinkowskiFunctional(1,1,0)] = 0.4
    D = MinkowskiDistribution(1, 1, 1, p_counter)

    @test 0.1 ≈ D.α[MinkowskiFunctional(0,0,0)]
    @test 0.3 ≈ D.α[MinkowskiFunctional(0,1,0)]
    @test 0.6 ≈ D.α[MinkowskiFunctional(1,0,0)]
    @test 1.0 ≈ D.α[MinkowskiFunctional(1,1,0)]

    D_A = marginalize(D, (:P, :χ))
    @test 0.3 ≈ D_A.p[(A =0,)]
    @test 0.7 ≈ D_A.p[(A =1,)]
    @test 0.3 ≈ D_A.α[(A =0,)]
    @test 1.0 ≈ D_A.α[(A =1,)]

    Ω = DensityOfStates(joinpath(SAMPLES_DIR, "structure_5x5"))
    D = MinkowskiDistribution(Ω, 10, 10)
    d_poisson = Distributions.Poisson(D.λ)
    p = 1 -cdf(d_poisson, D.ρ-1)
    D_A_binomial = Binomial(D.n^2, p)
    D_A = marginalize(D, (:P, :χ))

    ps_equal = Bool[]
    for x in 1:25
        push!(ps_equal, pdf(D_A_binomial, x) ≈ D_A.p[(A=x,)])
    end
    @test sum(ps_equal) == 25

    ps = pdf(D_A_binomial)
    αs_equal = Bool[]
    for x in 1:25
        push!(αs_equal, round(sum(ps[ps .<= pdf(D_A_binomial, x)]), digits=10) ≈ round(D_A.α[(A=x,)], digits=10))
    end
    @test sum(αs_equal) == 25
end