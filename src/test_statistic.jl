using StatsBase
using DataStructures, Distributions
using MinkowskiFunctionals

struct AreaDistributionX
    n::Int
    λ::Float64
    ρ::Int
    p::Binomial{Float64}
    pvalues::Dict{Int64, Float64}
end

function AreaDistributionX(n, λ, ρ)
    d_poisson = Distributions.Poisson(λ)
    p = 1 - cdf(d_poisson, ρ-1)
    d_A = Binomial(n, p)

    xs = support(d_A)
    ps = pdf(d_A)

    pvalues = [sum(ps[ps .<= p]) for p in ps]
    pvalues = Dict(x => p for (x,p) in zip(xs, pvalues))

    AreaDistributionX(n, λ, ρ, d_A, pvalues)
end

function AreaDistributionX(n, λ, ρ, pvalues)
    d_poisson = Distributions.Poisson(λ)
    p = 1 - cdf(d_poisson, ρ-1)
    d_A = Binomial(n, p)

    xs = support(d_A)
    ps = pdf(d_A)

    AreaDistributionX(n, λ, ρ, d_A, pvalues)
end

function compatibility(d::AreaDistributionX, f::Int64)
    return d.pvalues[f]
end

function compatibility(d::AreaDistributionX, x::CountsMap)
    bw_map = BWMap(x, d.ρ)
    return compatibility(d, sum(bw_map.pixels))
end

function append!(h5f::HDF5.File, x::AreaDistributionX)
    if !("thresholds" in keys(h5f))
        g_thresholds = create_group(h5f, "thresholds")
        g_thresholds["n"] = x.n
        g_thresholds["λ"] = x.λ
    end

    if !("$(x.ρ)" in keys(h5f["thresholds"]))
        g = create_group(h5f, "thresholds/$(x.ρ)")
        g["pvalues"] = [(A, p) for (A, p) in x.pvalues]
        g["p"] = x.p.p
    end
end

function AreaDistributionX(h5f::HDF5.File, ρ)
    n = read(h5f["thresholds/n"])
    λ = read(h5f["thresholds/λ"])
    p_black = read(h5f["thresholds/$(ρ)/p"])
    pvalues = read(h5f["thresholds/$(ρ)/pvalues"])
    d_pvalues = Dict(k => v for (k, v) in pvalues)
    p = Binomial(n, p_black)

    AreaDistributionX(n, λ, ρ, p, d_pvalues)
end

struct ECCDF
    λ::Float64
    L::Int64
    eccdf::Dict{Float64, Float64}
end

function ECCDF(λ::Float64, L::Int64, ecdf::T, xs) where {T <: ECDF}
    ys = 1 .- ecdf.(xs)

    ECCDF(λ, L, Dict(x=>y for (x,y) in zip(xs, ys)))
end

function (eccdf::ECCDF)(x::Float64)
    ks = sort(collect(keys(eccdf.eccdf)))
    closest_x = ks[searchsortednearest(ks, x)]
    return eccdf.eccdf[closest_x]
end

function compatibility(eccdf::ECCDF, dd::DefaultDict{Int64, AreaDistributionX, Int64}, x::Union{CountsMap, Matrix{Int64}})
    ts = calc_ts!(dd, x)

    return eccdf(ts)
end

function calc_ts!(dd::DefaultDict{Int64, AreaDistributionX, Int64}, x::Union{CountsMap, Matrix{Int64}})
    L, _ = size(x)
    ρs = get_thresholds(x)
    for ρ in ρs
        kys = keys(dd)
        if !(ρ in kys)
            dd[ρ] = AreaDistributionX(L^2, b, ρ)
        end
    end
    summed_ts = 0.0
    for ρ in ρs
        summed_ts += -log10(compatibility(dd[ρ], x))
    end

    return summed_ts
end

function sample_pvalues(dd, b, N)
    sum_pvalues = zeros(N)
    L = Int(sqrt(dd[1].n))
    @threads for i in 1:N
        sum_pvalues[i] = calc_ts!(dd, CountsMap(L, b))
    end

    return sum_pvalues
end

function sample(W, N, λ)
    d = Dict(1 => AreaDistributionX(W^2, λ, 1))
    dd = DefaultDict(0, d)
    ts = sample_pvalues(dd, λ, N)
    e_cdf = ecdf(ts)
    xs = 1:0.01:maximum(ts)
    ys = p2σ.(1 .- e_cdf.(xs))
    return xs, ys, dd
end
