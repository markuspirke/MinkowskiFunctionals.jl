using Test
using MinkowskiFunctionals
using Distributions

const SAMPLES_DIR = joinpath(@__DIR__, "samples")

@testset "minkowski_map" begin
    mink_map = MinkowskiMap(1.0, 2, 3, ones(3, 3))
    @test (3, 3) == size(mink_map)
    @test 1.0 == mink_map[1, 1]

    counts_map = CountsMap([1 1 1; 1 2 1; 1 1 1])
    L = 3
    b = Background(ones(3, 3))

    p_black1 = 1 - cdf(Poisson(1.0), 1-1) # prob of black pixel for rho = 1
    D1 = Binomial(L^2, p_black1) # binomial
    A1 = 9 # for treshold of 1
    p_A1 = pdf(D1, A1) # probability for A = 9
    pvalue1 = sum(pdf(D1)[pdf(D1) .<= p_A1]) # sum of probabilities less likeley than A=9
    s1 = 1.0

    p_black2 = 1 - cdf(Poisson(1.0), 2-1)
    D2 = Binomial(L^2, p_black2)
    A2 = 1 # for treshold of 1
    p_A2 = pdf(D2, A2) # probability for A = 9
    pvalue2 = sum(pdf(D2)[pdf(D2) .<= p_A2]) # sum of probabilities less likeley than A=9
    s2 = -1.0 # less pixels black than expected


    pvalue = minimum([pvalue1, pvalue2]) # 1 is smaller
    pvalue = 1 - (1 - pvalue)^2 #
    pvalue *= s1

    @test pvalue ≈ MinkowskiMap(counts_map, b, L).pixels[1, 1]
    @test pvalue ≈ MinkowskiMap(counts_map, 1.0, L).pixels[1, 1]

    Ω = DensityOfStates(3)
    D1 = MinkowskiDistribution(Ω, 1.0, 1)
    f1 = MinkowskiFunctional(9, 12, 1)
    p1 = pdf(D1, f1)
    pvalue1 = sum(pdf(D1)[pdf(D1) .<= p1])

    D2 = MinkowskiDistribution(Ω, 1.0, 2)
    f2 = MinkowskiFunctional(1, 4, 1)
    p2 = pdf(D2, f2)
    pvalue2 = sum(pdf(D2)[pdf(D2) .<= p2])
    pvalue = minimum([pvalue1, pvalue2]) # 1 is smaller
    pvalue = 1 - (1 - pvalue)^2 #
    pvalue *= s1

    @test pvalue ≈ MinkowskiMap(counts_map, b, Ω).pixels[1, 1]
    @test pvalue ≈ MinkowskiMap(counts_map, 1.0, Ω).pixels[1, 1]
end
