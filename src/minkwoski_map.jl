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
    αs = zeros(n - 2l, m - 2l)
    for j in l+1:m-l
        αs_ρ = zeros(length(ρs))
        for i in l+1:n-l
            for (k, ρ) in enumerate(ρs)
                αs_ρ[k] = compatibility(x[i-l:i+l, j-l:j+l], h0_distributions[k])
                # bw_map = BWMap(x[i-l:i+l, j-l:j+l], ρ)
                # functional = MinkowskiFunctional(bw_map.pixels)
                # deviation_strengths[k] = h0_distributions[k].α[functional]
            end
            @inbounds αs[i-l, j-l] = minimum(αs_ρ)
        end
    end

    return MinkowskiMap(λ, ρs, L, αs)
end

"""
    function compatibility(x::Matrix{Int64}, h0::MinkowskiDistribution)

Given some matrix of counts, this calculates the α-value for the center pixel
based on the all Minkowski functionals.
"""
function compatibility(x::Matrix{Int64}, h0::MinkowskiDistribution)
    bw_map = BWMap(x, h0.ρ)
    functional = MinkowskiFunctional(bw_map.pixels)

    h0.α[functional]
end

"""
    function MinkowskiMap(x::CountsMap, h0s::Vector{AreaDistribution})

This calculates a Minkowski Map only based on the area functional A.
This is possible for arbitrary window sizes, as the distribution of
the area functional is just a Binomial distribution.
"""
function MinkowskiMap(x::CountsMap, h0s::Vector{AreaDistribution})
    m, n = size(x.pixels)
    ρs = [s.ρ for s in h0s]
    λ = h0s[1].λ
    L = h0s[1].n
    l = floor(Int, L/2)
    αs = zeros(n - 2l, m - 2l)
    for j in l+1:m-l
        αs_ρ = zeros(length(ρs))
        for i in l+1:n-l
            for (k, ρ) in enumerate(ρs)
                αs_ρ[k] = compatibility(x[i-l:i+l, j-l:j+l], h0s[k])
            end
            @inbounds αs[i-l, j-l] = minimum(αs_ρ)
        end
    end

    return MinkowskiMap(λ, ρs, L, αs)
end

"""
    function compatibility(x::Matrix{Int64}, h0::AreaDistribution)

Given some matrix of counts, this calculates the α-value for the center pixel
based on the area functional.
"""
function compatibility(x::Matrix{Int64}, h0::AreaDistribution)
    bw_map = BWMap(x, h0.ρ)
    functional = MinkowskiFunctional(bw_map.pixels)

    h0.α[functional.A]
end

function MinkowskiMap(x::CountsMap, h0_distributions::Vector{MinkowskiDistribution}, fields)
    m, n = size(x.pixels)
    ρs = [s.ρ for s in h0_distributions]
    λ = h0_distributions[1].λ
    L = h0_distributions[1].n
    l = floor(Int, L/2)
    αs = zeros(n - 2l, m - 2l)
    for j in l+1:m-l
        αs_ρ = zeros(length(ρs))
        for i in l+1:n-l
            for (k, ρ) in enumerate(ρs)
                bw_map = BWMap(x[i-l:i+l, j-l:j+l], ρ)
                functional = MinkowskiFunctional(bw_map.pixels)
                reduced_functional = reduce_functional(functional, fields)
                αs_ρ[k] = h0_distributions[k].α[reduced_functional]
            end
            @inbounds αs[i-l, j-l] = minimum(αs_ρ)
        end
    end

    return MinkowskiMap(λ, ρs, L, αs)
end

"""
    function minkowski_map_A(x::CountsMap, λ, ρs, L)

This calculates a Minkowski Map only based on the area functional A.
This is possible for arbitrary window sizes, as the distribution of
the area functional is just a Binomial distribution.
"""
function minkowski_map_A(x::CountsMap, λ, ρs, L)
    m, n = size(x.pixels)
    dAs = [AreaDistribution(L, λ, ρ) for ρ in ρs]

    l = floor(Int, L/2)
    αs = zeros(n - 2l, m - 2l) # DS not appropirate name as this is not dev strength
    for j in l+1:m-l
        αs_ρ = zeros(length(ρs))
        for i in l+1:n-l
            for (k, ρ) in enumerate(ρs)
                bw_map = BWMap(x[i-l:i+l, j-l:j+l], ρ)
                functional = MinkowskiFunctional(bw_map.pixels)
                A = functional.A
                αs_ρ[k] = αs_ρ[k][A]
            end
            @inbounds αs[i-l, j-l] = minimum(αs_ρ)
        end
    end

    return MinkowskiMap(λ, ρs, L, αs)
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
    N = sum(round_mask)

    d_poisson = Distributions.Poisson(λ)
    ps = [1 - cdf(d_poisson, ρ-1) for ρ in ρs]
    ds_A = [Binomial(N, p) for p in ps]
    αs_A = [get_αs_binomial(d) for d in ds_A]

    αs = zeros(n - 2l, m - 2l) # DS not appropirate name as this is not dev strength
    for j in l+1:m-l
        αs_ρ = zeros(length(ρs))
        for i in l+1:n-l
            for (k, ρ) in enumerate(ρs)
                bw_map = BWMap(x[i-l:i+l, j-l:j+l] .* round_mask, ρ)
                functional = MinkowskiFunctional(bw_map.pixels)
                A = functional.A
                αs_ρ[k] = αs_A[k][A]
            end
            @inbounds αs[i-l, j-l] = minimum(αs_ρ)
        end
    end

    return MinkowskiMap(λ, ρs, L, αs)
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


function correct_trials(mink_map::MinkowskiMap, n_trials)
    pixels = mink_map.pixels * n_trials
    return  MinkowskiMap(mink_map.λ, mink_map.ρ, mink_map.L, pixels)
end

function p2σ(mink_map::MinkowskiMap)
    pixels = mink_map.pixels
    pixels[pixels .> 1] .= 1.0
    pixels = p2σ.(pixels)
    return  MinkowskiMap(mink_map.λ, mink_map.ρ, mink_map.L, pixels)
end
