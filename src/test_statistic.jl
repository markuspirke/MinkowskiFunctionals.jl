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
    y = eccdf.eccdf[closest_x]
    return y > 0.0 ? y : 1/eccdf.N
end

function compatibility(eccdf::ECCDF, dd::DefaultDict{Int64, T, Int64}, x::Union{CountsMap, Matrix{Int64}}) where {T<:AbstractMinkowskiDistribution}
    ts = calc_ts(dd, x)

    return eccdf(ts)
end

function compatibility(ecdf::T, dd::DefaultDict{Int64, S, Int64}, x::Union{CountsMap, Matrix{Int64}}) where {T <: ECDF, S <: AbstractMinkowskiDistribution}
    ts = calc_ts(dd, x)
    if 1.0 - ecdf(ts) > 0.0
        return 1.0 - ecdf(ts)
    else
        return length(e_cdf.sorted_values)
    end
end

function calc_ts(dd::DefaultDict{Int64, T, Int64}, x::Union{CountsMap, Matrix{Int64}}) where {T<:AbstractMinkowskiDistribution}
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

function MinkowskiMap(x::CountsMap, mink_ds::DefaultDict{Int64, T, Int64}, eccdf::ECCDF) where {T<:AbstractMinkowskiDistribution}
    λ = mink_ds[first(eachindex(mink_ds))].λ
    m, n = size(x)
    L = window_size(mink_ds[1])
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

function MinkowskiMap(x::CountsMap, mink_ds::DefaultDict{Int64, S, Int64}, eccdf::T) where {S<:AbstractMinkowskiDistribution, T<:ECDF}
    λ = mink_ds[first(eachindex(mink_ds))].λ
    m, n = size(x)
    L = window_size(mink_ds[1])
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
