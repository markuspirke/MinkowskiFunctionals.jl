"""
This generates a Minkowski map from a given CountsMap.
This also needs distributions of Minkowski functionals
and one needs to decide which functional is used.
"""
struct MinkowskiMap
    λ
    ρ
    L
    pixels
end

Base.size(x::MinkowskiMap) = size(x.pixels)
Base.getindex(x::MinkowskiMap, i, j) = x.pixels[i, j]

"""
    function MinkowskiMap(x, h0_distributions, F)

Given an counts map x and the distributions for different tresholds ρ
this calculates a Minkowski map for a given functional F.
"""
function MinkowskiMap(x, h0_distributions, F)
    m, n = size(x.pixels)
    ρs = [s.ρ for s in h0_distributions]
    λ = h0_distributions[1].λ
    L = h0_distributions[1].n
    l = floor(Int, L/2)
    Ds = zeros(n - 2l, m - 2l)
    for j in l+1:m-l
        for i in l+1:n-l
            deviation_strengths = zeros(length(ρs))
            for (k, ρ) in enumerate(ρs)
                bw_map = BWMap(x[i-l:i+l, j-l:j+l], ρ)
                functional = MinkowskiFunctional(bw_map.pixels)
                deviation_strengths[k] = deviation_strength(getfield(h0_distributions[k], F), getfield(functional, F))
            end
            Ds[i-l, j-l] = maximum(deviation_strengths)
        end
    end

    return MinkowskiMap(λ, ρs, L, Ds)
end

"""
    function MinkowskiMap(x::CountsMap, h0_distributions::Vector{MinkowskiDistribution})

Given an counts map x and the distributions for different tresholds ρ
this calculates a Minkowski map based on all three Minkowski functionals.
"""
function MinkowskiMap(x::CountsMap, h0_distributions::Vector{MinkowskiDistribution})
    m, n = size(x.pixels)
    ρs = [s.ρ for s in h0_distributions]
    λ = h0_distributions[1].λ
    L = h0_distributions[1].n
    l = floor(Int, L/2)
    Ds = zeros(n - 2l, m - 2l)
    for j in l+1:m-l
        deviation_strengths = zeros(length(ρs))
        for i in l+1:n-l
            for (k, ρ) in enumerate(ρs)
                bw_map = BWMap(x[i-l:i+l, j-l:j+l], ρ)
                functional = MinkowskiFunctional(bw_map.pixels)
                deviation_strengths[k] = h0_distributions[k].α[functional]
            end
            @inbounds Ds[i-l, j-l] = minimum(deviation_strengths)
        end
    end

    return MinkowskiMap(λ, ρs, L, Ds)
end

function MinkowskiMap(x::CountsMap, h0_distributions::Vector{MinkowskiDistribution}, fields)
    m, n = size(x.pixels)
    ρs = [s.ρ for s in h0_distributions]
    λ = h0_distributions[1].λ
    L = h0_distributions[1].n
    l = floor(Int, L/2)
    Ds = zeros(n - 2l, m - 2l)
    for j in l+1:m-l
        deviation_strengths = zeros(length(ρs))
        for i in l+1:n-l
            for (k, ρ) in enumerate(ρs)
                bw_map = BWMap(x[i-l:i+l, j-l:j+l], ρ)
                functional = MinkowskiFunctional(bw_map.pixels)
                reduced_functional = reduce_functional(functional, fields)
                deviation_strengths[k] = h0_distributions[k].α[reduced_functional]
            end
            @inbounds Ds[i-l, j-l] = minimum(deviation_strengths)
        end
    end

    return MinkowskiMap(λ, ρs, L, Ds)
end

"""
    function minkowski_map_A(x::CountsMap, λ, ρs, L)

This calculates a Minkowski Map only based on the area functional A.
This is possible for arbitrary window sizes, as the distribution of
the area functional is just a Binomial distribution.
"""
function minkowski_map_A(x::CountsMap, λ, ρs, L)
    m, n = size(x.pixels)
    d_poisson = Distributions.Poisson(λ)
    ps = [1 - cdf(d_poisson, ρ-1) for ρ in ρs]
    ds_A = [Binomial(L^2, p) for p in ps]
    αs_A = [get_αs_binomial(d) for d in ds_A]

    l = floor(Int, L/2)
    Ds = zeros(n - 2l, m - 2l) # DS not appropirate name as this is not dev strength
    for j in l+1:m-l
        αs_ρ = zeros(length(ρs))
        for i in l+1:n-l
            for (k, ρ) in enumerate(ρs)
                bw_map = BWMap(x[i-l:i+l, j-l:j+l], ρ)
                functional = MinkowskiFunctional(bw_map.pixels)
                A = functional.A
                αs_ρ[k] = αs_A[k][A]
            end
            @inbounds Ds[i-l, j-l] = minimum(αs_ρ)
        end
    end

    return MinkowskiMap(λ, ρs, L, Ds)
end

"""
    function get_αs_binomial(d::Binomial{Float64})

A small helper function to precalculate a lookup table for
hypothesis testing based on a Binomial distribution.
"""
function get_αs_binomial(d::Binomial{Float64})
    ps = pdf(d)
    αs = Dict(i-1 => sum(ps[ps[i] .>= ps]) for i in 1:length(ps))

    return αs
end

"""
    function minkowski_map_A_round(x::CountsMap, λ, ρs, L)

This calculates a Minkowski Map only based on the area functional A with
a circular kernel.
This is possible for arbitrary window sizes, as the distribution of
the area functional is just a Binomial distribution.
"""
function minkowski_map_A_round(x::CountsMap, λ, ρs, L)
    m, n = size(x.pixels)
    r = floor(Int, L/2)
    l = floor(Int, L/2)
    round_mask = [sqrt((i - r - 1)^2 + (j - r - 1)^2) <= r for i in 1:L, j in 1:L]
    @show round_mask
    N = sum(round_mask)
    @show N

    d_poisson = Distributions.Poisson(λ)
    ps = [1 - cdf(d_poisson, ρ-1) for ρ in ρs]
    ds_A = [Binomial(N, p) for p in ps]
    αs_A = [get_αs_binomial(d) for d in ds_A]

    Ds = zeros(n - 2l, m - 2l) # DS not appropirate name as this is not dev strength
    for j in l+1:m-l
        αs_ρ = zeros(length(ρs))
        for i in l+1:n-l
            for (k, ρ) in enumerate(ρs)
                bw_map = BWMap(x[i-l:i+l, j-l:j+l] .* round_mask, ρ)
                functional = MinkowskiFunctional(bw_map.pixels)
                A = functional.A
                αs_ρ[k] = αs_A[k][A]
            end
            @inbounds Ds[i-l, j-l] = minimum(αs_ρ)
        end
    end

    return MinkowskiMap(λ, ρs, L, Ds)
end

"""
    function deviation_strength(h0_distribution::T, x) where {T<:DiscreteDistribution}

This takes a discrete distribution and calulates the deviation strength for
the event X, where X needs to be in the support (sample space) of the distribution.
"""
function deviation_strength(h0_distribution::T, x) where {T<:DiscreteDistribution}
    p = pdf(h0_distribution, x)
    mask = probs(h0_distribution) .<= p
    return -log10(sum(probs(h0_distribution)[mask]))
end
