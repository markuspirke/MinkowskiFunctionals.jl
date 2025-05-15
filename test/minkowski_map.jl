using Test
using MinkowskiFunctionals
using Distributions

const SAMPLES_DIR = joinpath(@__DIR__, "samples")

@testset "minkowski_map" begin
    mink_map = MinkowskiMap(ones(3, 3))
    @test (3, 3) == size(mink_map)
    @test 1.0 == mink_map[1, 1]

    a = 1e-20
    @test 2*a ≈ correct_trials(a, 2)
    @test 3*a ≈ correct_trials(a, 3)
    b = 0.1
    @test 1 - (1 - b)^2 ≈ correct_trials(b, 2)
    @test 1 - (1 - b)^3 ≈ correct_trials(b, 3)


    counts_map = CountsMap([1 1 1; 2 3 1; 1 1 1])
    L = 3
    b = Background(ones(3, 3))

    p_black1 = 1 - cdf(Poisson(1.0), 2-1) # prob of black pixel for rho = 1
    D1 = Binomial(L^2, p_black1) # binomial
    A1 = 2 # for treshold of 1
    p_A1 = pdf(D1, A1) # probability for A = 9
    pvalue1 = sum(pdf(D1)[pdf(D1) .<= p_A1]) # sum of probabilities less likeley than A=9
    s1 = -1.0

    p_black2 = 1 - cdf(Poisson(1.0), 3-1)
    D2 = Binomial(L^2, p_black2)
    A2 = 1 # for treshold of 1
    p_A2 = pdf(D2, A2) # probability for A = 9
    pvalue2 = sum(pdf(D2)[pdf(D2) .<= p_A2]) # sum of probabilities less likeley than A=9
    s2 = 1.0 # less pixels black than expected

    pvalues = [pvalue1, pvalue2]
    signs = [s1, s2]
    idx = argmin(pvalues) # 1 is smaller
    pvalue = 1 - (1 - pvalues[idx])^2 #
    pvalue *= signs[idx]

    @test pvalue ≈ MinkowskiMap(counts_map, b, L).pixels[1, 1]
    @test pvalue ≈ MinkowskiMap(counts_map, 1.0, L).pixels[1, 1]
    mask = ones(3, 3) .== 1.0
    @test pvalue ≈ MinkowskiMap(counts_map, b, mask).pixels[1, 1]

    r = floor(Int64, L/2)
    mask = [sqrt((i - r - 1)^2 + (j - r - 1)^2) <= r for i in 1:L, j in 1:L]

    p_black1 = 1 - cdf(Poisson(1.0), 2-1) # prob of black pixel for rho = 1
    D1 = Binomial(5, p_black1) # binomial
    A1 = 2 # for treshold of 1
    p_A1 = pdf(D1, A1) # probability for A = 9
    pvalue1 = sum(pdf(D1)[pdf(D1) .<= p_A1]) # sum of probabilities less likeley than A=9
    s1 = 1.0

    p_black2 = 1 - cdf(Poisson(1.0), 3-1)
    D2 = Binomial(5, p_black2)
    A2 = 1 # for treshold of 1
    p_A2 = pdf(D2, A2) # probability for A = 9
    pvalue2 = sum(pdf(D2)[pdf(D2) .<= p_A2]) # sum of probabilities less likeley than A=9
    s2 = 1.0 # less pixels black than expected

    pvalues = [pvalue1, pvalue2]
    signs = [s1, s2]
    idx = argmin(pvalues) # 1 is smaller
    pvalue = 1 - (1 - pvalues[idx])^2 #
    pvalue *= signs[idx]
    # pvalue = minimum([pvalue1, pvalue2]) # 1 is smaller
    # pvalue = 1 - (1 - pvalue)^2 #
    # pvalue *= s1
    @test pvalue ≈ MinkowskiMap(counts_map, b, mask).pixels[1, 1]

    Ω = DensityOfStates(3)
    D1 = MinkowskiDistribution(Ω, 1.0, 2)
    f1 = MinkowskiFunctional(2, 6, 1)
    p1 = pdf(D1, f1)
    pvalue1 = sum(pdf(D1)[pdf(D1) .<= p1])

    D2 = MinkowskiDistribution(Ω, 1.0, 3)
    f2 = MinkowskiFunctional(1, 4, 1)
    p2 = pdf(D2, f2)
    pvalue2 = sum(pdf(D2)[pdf(D2) .<= p2])


    pvalues = [pvalue1, pvalue2]
    signs = [s1, s2]
    idx = argmin(pvalues) # 1 is smaller
    pvalue = 1 - (1 - pvalues[idx])^2 #
    pvalue *= signs[idx]

    @test pvalue ≈ MinkowskiMap(counts_map, b, Ω).pixels[1, 1]
    @test pvalue ≈ MinkowskiMap(counts_map, 1.0, Ω).pixels[1, 1]

    counts = zeros(3, 3)
    counts_map = CountsMap(counts)
    background = Background(zeros(3, 3))
    mink_map = MinkowskiMap(counts_map, background, Ω)
    @test 1.0 ≈ mink_map[1, 1]

    background = Background([1.0 1.0 100000.0; 1.0 1.0 9.0; 1.0 1.0 9.0])
    counts_map = CountsMap(3, 2.0)
    mink_map = MinkowskiMap(counts_map, background, Ω)


    Ω = DensityOfStates(3)
    counts_map = CountsMap(3, 2.0)
    background = Background(2*ones(3, 3))
    background_float = 2.0
    @test MinkowskiMap(counts_map, background, Ω).pixels ≈ MinkowskiMap(counts_map, background_float, Ω).pixels

    ρs = get_thresholds(counts_map)
    mink_ds = Dict(ρ => MinkowskiDistribution(Ω, background_float, ρ) for ρ in ρs)
    @test MinkowskiMap(counts_map, background, Ω).pixels ≈ MinkowskiMap(counts_map, background_float, mink_ds).pixels

    mink_ds = Dict(ρ => AreaDistribution(3^2, background_float, ρ) for ρ in ρs)
    @test MinkowskiMap(counts_map, background, 3).pixels ≈ MinkowskiMap(counts_map, background_float, mink_ds).pixels

    # check if correction works at boundaries
    as = Float64[]
    bs = Float64[]
    for _ in 1:100
        background = zeros(32, 32)
        background[:, 1:16] .= 2.0
        background = Background(background)
        counts_map = CountsMap(background)
        mink_map = MinkowskiMap(counts_map, background, 9)
        x = p2σ(mink_map)
        push!(as, mean(x[:, 1:6]))
        push!(bs, mean(x[:, 7:12]))
    end

    @test 0.1 > abs(mean(as) - mean(bs))

    Ω = DensityOfStates(joinpath(SAMPLES_DIR, "structure_5x5"))
    as = Float64[]
    bs = Float64[]
    for _ in 1:50
        background = zeros(12, 12)
        background[:, 1:6] .= 2.0
        background = Background(background)
        counts_map = CountsMap(background)
        mink_map = MinkowskiMap(counts_map, background, Ω)
        x = p2σ(mink_map)
        push!(as, mean(x[:, 1:2]))
        push!(bs, mean(x[:, 3:4]))
    end

    @test 0.3 > abs(mean(as) - mean(bs))

    as = Float64[]
    bs = Float64[]
    for _ in 1:50
        background = ones(12, 12)
        background[:, 1:6] .= 3.0
        background = Background(background)
        counts_map = CountsMap(background)
        mink_map = MinkowskiMap(counts_map, background, Ω)
        x = p2σ(mink_map)
        push!(as, mean(x[:, 5:6]))
        push!(bs, mean(x[:, 7:8]))
    end

    @test 0.3 > abs(mean(as) - mean(bs))

    # SAME TEST FOR AREA ONLY
    as = Float64[]
    bs = Float64[]
    for _ in 1:100
        background = zeros(12, 12)
        background[:, 1:6] .= 2.0
        background = Background(background)
        counts_map = CountsMap(background)
        mink_map = MinkowskiMap(counts_map, background, 5)
        x = p2σ(mink_map)
        push!(as, mean(x[:, 1:2]))
        push!(bs, mean(x[:, 3:4]))
    end

    @test 0.1 > abs(mean(as) - mean(bs))

    as = Float64[]
    bs = Float64[]
    for _ in 1:100
        background = ones(12, 12)
        background[:, 1:6] .= 3.0
        background = Background(background)
        counts_map = CountsMap(background)
        mink_map = MinkowskiMap(counts_map, background, 5)
        x = p2σ(mink_map)
        push!(as, mean(x[:, 5:6]))
        push!(bs, mean(x[:, 7:8]))
    end

    @test 0.1 > abs(mean(as) - mean(bs))
end
