using StatsBase
using DataStructures, Distributions
using MinkowskiFunctionals

struct ECCDF
    λ::Float64
    L::Int64
    N::Int64
    eccdf::SortedDict{Float64, Float64}
end

function ECCDF(λ::Float64, L::Int64, ecdf::T, N) where {T <: ECDF}
    xs = vcat(range(minimum(ecdf), maximum(ecdf), N))
    xs[end] -= 1e-8
    ys = 1 .- ecdf.(xs)

    ECCDF(λ, L, length(ecdf.sorted_values), SortedDict(Dict(x=>y for (x,y) in zip(xs, ys))))
end

function ECCDF(λ::Float64, L::Int64, N, n)
    d = Dict(1 => AreaDistribution(L^2, λ, 1))
    dd = DefaultDict(0, d)
    ts = [calc_ts!(dd, CountsMap(L, λ)) for _ in 1:N]
    e_cdf = ecdf(ts)
    eccdf = ECCDF(λ, L, e_cdf, n)

    return eccdf
end

function (eccdf::ECCDF)(x::Float64)
    ks = collect(keys(eccdf.eccdf))
    idx = searchsortedlast(ks, x)
    idx == 0 && return 1.0
    closest_x = ks[idx]

    return eccdf.eccdf[closest_x]
end

function compatibility(eccdf::ECCDF, dd::DefaultDict{Int64, AreaDistribution, Int64}, x::Union{CountsMap, Matrix{Int64}})
    ts = calc_ts(dd, x)

    return eccdf(ts)
end

function compatibility(eccdf::ECCDF, dd::DefaultDict{Int64, MinkowskiDistribution, Int64}, x::Union{CountsMap, Matrix{Int64}})
    ts = calc_ts(dd, x)

    return eccdf(ts)
end

function compatibility(ecdf::T, dd::DefaultDict{Int64, MinkowskiDistribution, Int64}, x::Union{CountsMap, Matrix{Int64}}) where {T <: ECDF}
    ts = calc_ts(dd, x)
    return 1.0 - ecdf(ts)
end

function calc_ts(dd::DefaultDict{Int64, MinkowskiDistribution, Int64}, x::Union{CountsMap, Matrix{Int64}})
    ρs = get_thresholds(x)
    summed_ts = 0.0
    for ρ in ρs
        summed_ts += -log10(compatibility(dd[ρ], x))
    end

    return summed_ts
end

function calc_ts!(dd::DefaultDict{Int64, AreaDistribution, Int64}, x::Union{CountsMap, Matrix{Int64}})
    b = dd[first(eachindex(dd))].λ
    L, _ = size(x)
    ρs = get_thresholds(x)
    for ρ in ρs
        kys = keys(dd)
        if !(ρ in kys)
            dd[ρ] = AreaDistribution(L^2, b, ρ)
        end
    end
    summed_ts = 0.0
    for ρ in ρs
        summed_ts += -log10(compatibility(dd[ρ], x))
    end

    return summed_ts
end

function calc_ts(dd::DefaultDict{Int64, AreaDistribution, Int64}, x::Union{CountsMap, Matrix{Int64}})
    ρs = get_thresholds(x)
    summed_ts = 0.0
    for ρ in ρs
        summed_ts += -log10(compatibility(dd[ρ], x))
    end

    return summed_ts
end

function calc_ts(dd::DefaultDict{Int64, MinkowskiDistribution, Int64}, x::Union{CountsMap, Matrix{Int64}})
    ρs = get_thresholds(x)
    summed_ts = 0.0
    for ρ in ρs
        summed_ts += -log10(compatibility(dd[ρ], x))
    end

    return summed_ts
end

function update_distributions!(dd::DefaultDict{Int64, AreaDistribution, Int64}, x::CountsMap)
    b = dd[first(eachindex(dd))].λ
    n = dd[first(eachindex(dd))].n
    ρs = get_thresholds(x)
    for ρ in ρs
        kys = keys(dd)
        if !(ρ in kys)
            dd[ρ] = AreaDistribution(n, b, ρ)
        end
    end
end

function MinkowskiMap(x::CountsMap, mink_ds::DefaultDict{Int64, AreaDistribution, Int64}, eccdf::ECCDF)
    λ = mink_ds[first(eachindex(mink_ds))].λ
    m, n = size(x)
    L = Int(sqrt(mink_ds[1].n))
    l = floor(Int, L/2)
    αs = zeros(n - 2l, m - 2l)
    signs = zeros(n - 2l, m - 2l)
    Threads.@threads for j in l+1:m-l
        for i in l+1:n-l
            local_counts = x[i-l:i+l, j-l:j+l]
            pvalue = compatibility(eccdf, mink_ds, CountsMap(local_counts))
            αs[i-l, j-l] = pvalue
            signs[i-l, j-l] = mean(local_counts) > λ ? 1.0 : -1.0
        end
    end
    MinkowskiMap(αs .* signs)
end

function MinkowskiMap(x::CountsMap, mink_ds::DefaultDict{Int64, MinkowskiDistribution, Int64}, eccdf::ECCDF)
    λ = mink_ds[first(eachindex(mink_ds))].λ
    m, n = size(x)
    L = mink_ds[1].n
    l = floor(Int, L/2)
    αs = zeros(n - 2l, m - 2l)
    signs = zeros(n - 2l, m - 2l)
    Threads.@threads for j in l+1:m-l
        for i in l+1:n-l
            local_counts = x[i-l:i+l, j-l:j+l]
            pvalue = compatibility(eccdf, mink_ds, CountsMap(local_counts))
            αs[i-l, j-l] = pvalue
            signs[i-l, j-l] = mean(local_counts) > λ ? 1.0 : -1.0
        end
    end
    MinkowskiMap(αs .* signs)
end

function MinkowskiMap(x::CountsMap, mink_ds, eccdf::T) where {T <: ECDF}
    λ = mink_ds[first(eachindex(mink_ds))].λ
    m, n = size(x)
    L = mink_ds[1].n
    l = floor(Int, L/2)
    αs = zeros(n - 2l, m - 2l)
    signs = zeros(n - 2l, m - 2l)
    Threads.@threads for j in l+1:m-l
        for i in l+1:n-l
            local_counts = x[i-l:i+l, j-l:j+l]
            pvalue = compatibility(eccdf, mink_ds, CountsMap(local_counts))
            αs[i-l, j-l] = pvalue
            signs[i-l, j-l] = mean(local_counts) > λ ? 1.0 : -1.0
        end
    end
    MinkowskiMap(αs .* signs)
end
