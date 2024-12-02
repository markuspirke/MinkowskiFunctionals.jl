using MinkowskiFunctionals
using Distributions
using Test

@testset "distribution" begin

    n, λ, ρ = 3, 10, 10
    Ω = DensityOfStates(n)
    @test 2^(n^2) == sum(values(Ω.data))
    P = MinkowskiDistribution(Ω, λ, ρ)
    @test 1.0 ≈ sum(pdf(P))

    P_A = marginalize(P, :A)
    @test 1.0 ≈ sum(pdf(P_A))
    P_P = marginalize(P, :P)
    @test 1.0 ≈ sum(pdf(P_P))
    P_χ = marginalize(P, :χ)
    @test 1.0 ≈ sum(pdf(P_χ))

    p = 1 - cdf(Distributions.Poisson(λ), ρ-1)
    P_A_from_binomial = Binomial(n^2, p)

    @test support(P_A_from_binomial) == support(P_A)
    @test (pdf(P_A_from_binomial) .≈ pdf(P_A)) == ones(Int, 10)

    # MINKOWSKI DEVIATION STRENGTH
    p_counter = Accumulator{MinkowskiFunctional, Float64}()
    p_counter[MinkowskiFunctional(0,0,0)] = 0.15865525393145696
    p_counter[MinkowskiFunctional(1,0,0)] = 0.15865525393145696
    p_counter[MinkowskiFunctional(2,0,0)] = 0.6826894921370861

    D = MinkowskiDistribution(1, 1, 1, p_counter)
    @test compatibility(D, MinkowskiFunctional(0,0,0)) ≈ 1.0 # in sigma
    @test compatibility(D, MinkowskiFunctional(1,0,0)) ≈ 1.0 # in sigma
    @test compatibility(D, MinkowskiFunctional(2,0,0)) ≈ 0.0 # in sigma



    h5open(joinpath(SAMPLES_DIR, "test_lookuptable.h5"), "w") do h5f
        append!(h5f, D)
    end

    D_loaded = MinkowskiDistribution(joinpath(SAMPLES_DIR, "test_lookuptable.h5"), 1, 1, 1)
    @test D.n == D_loaded.n
    @test D.λ == D_loaded.λ
    @test D.ρ == D_loaded.ρ
    @test D.P == D_loaded.P
    @test D.σ == D_loaded.σ

end
