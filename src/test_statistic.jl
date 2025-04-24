using Base.Threads
using StatsBase
using DataStructures, Distributions
using MinkowskiFunctionals
import MinkowskiFunctionals: compatibility

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

# d = Dict(x => AreaDistributionX(L^2, b, x) for x in 1:Int(4*b))

function sum_over_thresholds!(dd::DefaultDict{Int64, AreaDistributionX, Int64}, x::CountsMap, b::Float64)
    L, _ = size(x)
    ρs = get_tresholds(x)
    for ρ in ρs
        kys = keys(dd)
        if !(ρ in kys)
            dd[ρ] = AreaDistributionX(L^2, b, ρ)
        end
    end
    sum_pvalue = 0.0
    for ρ in ρs
        sum_pvalue += -log10(compatibility(dd[ρ], x))
    end

    return sum_pvalue
end


function foo(dd, b, N)
    sum_pvalues = zeros(N)
    L = Int(sqrt(d[1].n))
    @threads for i in 1:N
        sum_pvalues[i] = sum_over_thresholds!(dd, CountsMap(L, b), b)
    end

    return sum_pvalues
end

# STORE ECDF AND AREA DISTRIBUTIONS
λs = [9.0, 10.0, 11.0]
W = 5
N = 1_000
d = Dict(1 => AreaDistributionX(W^2, λs[1], 1))
dd = DefaultDict(0, d)
h5f = h5open("foo.h5", "w")
g_ecdf = create_group(h5f, "ecdf")
g_pdf = create_group(h5f, "ecdf")
xs = 1:0.1:60
for λ in λs
    d = Dict(1 => AreaDistributionX(W^2, λ, 1))
    dd = DefaultDict(0, d)
    pvalues = foo(dd, λ, N)
    e_cdf = ecdf(pvalues)
    ys = p2σ.(1 .- e_cdf.(xs))
    g_ecdf["λ=$(λ)"] = [(x, y) for (x, y) in zip(xs, ys)]
    # create_dataset(g_ecdf, "λ=$(λ)", [(x, y) for (x, y) in zip(xs, ys)])
end
close(h5f)
h5f = h5open("foo.h5", "r")






# MINK MAPS WITH NEW TEST STATISTIC
λ = 10.0
L = 32
countsmap = CountsMap(L, λ)
d_background = Normal(λ, 0.1)
background = Background(rand(d_background, L, L))



function minkmap_new_ts(x::CountsMap, b::Background, L::Int64)
    m, n = size(x.pixels)
    l = floor(Int, L/2)
    αs = ones(n - 2l, m - 2l)
    signs = ones(n - 2l, m - 2l)
    λs = unique(round.(background.pixels, digits=2))


function MinkowskiMap(x::CountsMap, b::Background, L::Int64)
    m, n = size(x.pixels)
    l = floor(Int, L/2)
    αs = ones(n - 2l, m - 2l)
    signs = ones(n - 2l, m - 2l)
    for j in l+1:m-l
        for i in l+1:n-l
            if b.pixels[i, j] == 0.0
                continue
            end
            local_counts = x[i-l:i+l, j-l:j+l]
            local_background = b[i-l:i+l, j-l:j+l]
            correction!(local_counts, local_background, b[i, j])
            ρs = get_tresholds(local_counts)
            l_ρ = length(ρs)
            αs_ρ = ones(l_ρ)
            signs_ρ = zeros(l_ρ)
            for k in 1:length(ρs)
                mink_distribution = AreaDistribution(L^2, b.pixels[i, j], ρs[k])
                @inbounds αs_ρ[k] = compatibility(mink_distribution, local_counts)
                @inbounds signs_ρ[k] = get_sign(mink_distribution, local_counts)
            end
            idx = argmin(αs_ρ)
            α = αs_ρ[idx]
            @inbounds αs[i-l, j-l] = correct_trials(α, l_ρ)
            @inbounds signs[i-l, j-l] = signs_ρ[idx]
        end
    end
    return MinkowskiMap(αs .* signs)
end
