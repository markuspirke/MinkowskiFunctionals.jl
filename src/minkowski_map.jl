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

function p2σ(x::MinkowskiMap)
    p2σ.(abs.(x.pixels)) .* sign.(x.pixels)
end

Base.size(x::MinkowskiMap) = size(x.pixels)
Base.getindex(x::MinkowskiMap, i, j) = x.pixels[i, j]

function MinkowskiMap(x::CountsMap, b::Background, L)

    m, n = size(x)
    ρs = ceil(Int, minimum(b.pixels)):maximum(x.pixels)
    l_ρ = length(ρs)
    l = floor(Int, L/2)
    αs = zeros(n - 2l, m - 2l)
    signs = zeros(n - 2l, m - 2l)

    for j in l+1:m-l
        for i in l+1:n-l
            αs_ρ = zeros(l_ρ)
            signs_ρ = zeros(l_ρ)
            for (k, ρ) in enumerate(ρs)
                mink_distribution = AreaDistribution(L^2, b.pixels[i, j], ρ)
                αs_ρ[k] = compatibility(mink_distribution, x[i-l:i+l, j-l:j+l])
                signs_ρ[k] = get_sign(mink_distribution, x[i-l:i+l, j-l:j+l])
            end
            idx = argmin(αs_ρ)
            α = αs_ρ[idx]
            @inbounds αs[i-l, j-l] = 1 - (1 - α)^l_ρ
            @inbounds signs[i-l, j-l] = signs_ρ[idx]
        end
    end

    return MinkowskiMap(b, ρs, L, αs .* signs)
end


function MinkowskiMap(x::CountsMap, b::Float64, L)

    m, n = size(x.pixels)
    ρs = b:maximum(x.pixels)
    l_ρ = length(ρs)
    l = floor(Int, L/2)
    αs = zeros(n - 2l, m - 2l)
    signs = zeros(n - 2l, m - 2l)
    mink_ds = [AreaDistribution(L^2, b, ρ) for ρ in ρs]
    for j in l+1:m-l
        for i in l+1:n-l
            αs_ρ = zeros(l_ρ)
            signs_ρ = zeros(l_ρ)
            for (k, ρ) in enumerate(ρs)
                αs_ρ[k] = compatibility(mink_ds[k], x[i-l:i+l, j-l:j+l])
                signs_ρ[k] = get_sign(mink_ds[k], x[i-l:i+l, j-l:j+l])
            end
            idx = argmin(αs_ρ)
            α = αs_ρ[idx]
            @inbounds αs[i-l, j-l] = 1 - (1 - α)^l_ρ
            @inbounds signs[i-l, j-l] = signs_ρ[idx]
        end
    end

    return MinkowskiMap(b, ρs, L, αs .* signs)
end

get_sign(d::AreaDistribution, x::CountsMap) = d.p.p*d.p.n > sum(BWMap(x, d.ρ).pixels) ? -1.0 : 1.0
get_sign(d::AreaDistribution, x::Matrix{Int64}) = d.p.p*d.p.n > sum(BWMap(x, d.ρ).pixels) ? -1.0 : 1.0


function MinkowskiMap(x::CountsMap, b::Background, Ω::DensityOfStates)
    m, n = size(x.pixels)
    ρs = ceil(Int, minimum(b.pixels)):maximum(x.pixels)
    l_ρ = length(ρs)
    L = Ω.n
    l = floor(Int, L/2)
    αs = zeros(n - 2l, m - 2l)
    signs = zeros(n - 2l, m - 2l)
    for j in l+1:m-l
        for i in l+1:n-l
            αs_ρ = zeros(l_ρ)
            signs_ρ = zeros(l_ρ)
            for (k, ρ) in enumerate(ρs)
                mink_distribution = MinkowskiDistribution(Ω, b.pixels[i, j], ρ)
                αs_ρ[k] = compatibility(mink_distribution, x[i-l:i+l, j-l:j+l])
                signs_ρ[k] = get_sign(mink_distribution, x[i-l:i+l, j-l:j+l])
            end
            idx = argmin(αs_ρ)
            α = αs_ρ[idx]
            @inbounds αs[i-l, j-l] = 1 - (1 - α)^l_ρ
            @inbounds signs[i-l, j-l] = signs_ρ[idx]
        end
    end

    return MinkowskiMap(b, ρs, L, αs .* signs)
end


function MinkowskiMap(x::CountsMap, b::Float64, Ω::DensityOfStates)
    m, n = size(x.pixels)

    ρs = b:maximum(x.pixels)
    l_ρ = length(ρs)
    L = Ω.n
    l = floor(Int, L/2)
    αs = zeros(n - 2l, m - 2l)
    signs = zeros(n - 2l, m - 2l)
    mink_ds = [MinkowskiDistribution(Ω, b, ρ) for ρ in ρs]
    for j in l+1:m-l
        for i in l+1:n-l
            αs_ρ = zeros(l_ρ)
            signs_ρ = zeros(l_ρ)
            for (k, ρ) in enumerate(ρs)
                αs_ρ[k] = compatibility(mink_ds[k], x[i-l:i+l, j-l:j+l])
                signs_ρ[k] = get_sign(mink_ds[k], x[i-l:i+l, j-l:j+l])
            end
            idx = argmin(αs_ρ)
            α = αs_ρ[idx]
            @inbounds αs[i-l, j-l] = 1 - (1 - α)^l_ρ
            @inbounds signs[i-l, j-l] = signs_ρ[idx]
        end
    end

    return MinkowskiMap(b, ρs, L, αs .* signs)
end

get_sign(d::MinkowskiDistribution, x::CountsMap) = d.p_black*d.n^2 > sum(BWMap(x, d.ρ).pixels) ? -1.0 : 1.0
get_sign(d::MinkowskiDistribution, x::Matrix{Int64}) = d.p_black*d.n^2 > sum(BWMap(x, d.ρ).pixels) ? -1.0 : 1.0


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
