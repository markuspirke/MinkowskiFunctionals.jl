"""
    struct ECCDF

Datatype which stores the empirical cumulative distibution function
for summed p-values at different thresholds.
"""
struct ECCDF{T<:Union{Float64, Background}}
    λ::T
    L::Int64
    N::Int64
    F::ECDF
    # ts::Vector{Float64}
    # pvalues::Vector{Float64}
end

# """
#     function ECCDF(λ::Float64, L::Int64, ecdf::T, N) where {T <: ECDF}

# Given an already calculated ecdf (from stats base) return the ECCDF.
# """
# function ECCDF(λ::Float64, L::Int64, ecdf::T, N) where {T <: ECDF}
#     xs = vcat(range(minimum(ecdf), maximum(ecdf), N))
#     ys = 1 .- ecdf.(xs)

#     ECCDF(λ, L, length(ecdf.sorted_values), xs, ys)
# end

# """
#     function ECCDF(λ::Background, L::Int64, ecdf::T, N) where {T <: ECDF}

# Given an already calculated ecdf (from stats base) return the ECCDF.
# """
# function ECCDF(λ::Background, L::Int64, ecdf::T, N) where {T <: ECDF}
#     xs = vcat(range(minimum(ecdf), maximum(ecdf), N))
#     ys = 1 .- ecdf.(xs)

#     ECCDF(λ, L, length(ecdf.sorted_values), xs, ys)
# end

"""
    function (eccdf::ECCDF)(x::Float64)

Finds the nearest pvalue within the ecdf. A self implemented
binary search algorithm is used, which converges quickly.
"""
function (eccdf::ECCDF)(x::Float64)
    # idx = binary_search(eccdf.ts, x)
    # return idx < length(eccdf.ts) ? eccdf.pvalues[idx] : 1/eccdf.N
    pvalue = 1.0 - eccdf.F(x)
    return pvalue > 0.0 ? pvalue : 1/eccdf.N
end


"""
    function find_max_threshold(λ, L)

Find the maximum threshold that is needed to sufficiently calculate the ecdf for
an background with mean λ.
"""
function find_max_threshold(λ, L)
    ρ = 1
    p = 1.0
    while p > 1e-12 / L^2
        p, _ = gamma_inc(ρ, λ)
        ρ += 1
    end

    return ρ
end


"""
    function ECCDF(λ::Float64, L::Int64, N, n)

Calculates the ecdf for a system size of L at an average background λ with the Area functional.
"""
function ECCDF(λ::Float64, L::Int64, N)
    ρmax = find_max_threshold(λ, L)
    d = Dict(ρ => AreaDistribution(L^2, λ, ρ) for ρ in 1:ρmax)
    ts = zeros(N)
    Threads.@threads for i in 1:N
        ts[i] = calc_ts(d, CountsMap(L, λ))
    end
    e_cdf = ecdf(ts)
    return ECCDF(λ, L, N, e_cdf)
end

"""
    function calc_ts(lut::MinkowskiPValueLookup, x::Union{CountsMap, Matrix{Int64}}, λ::Float64)

Calculates the test statistic for a given counts map or matrix using a precomputed
`MinkowskiPValueLookup`. Replaces building a `Dict` of `MinkowskiDistribution`s per λ.
"""
function calc_ts(lut::MinkowskiPValueLookup, x::Union{CountsMap, Matrix{Int64}}, λ::Float64)
    pixels = x isa CountsMap ? x.pixels : x
    ρs = get_thresholds(pixels)
    summed_ts = 0.0
    for ρ in ρs
        p, _ = gamma_inc(ρ, λ)
        bw = BWMap(pixels, ρ)
        f  = MinkowskiFunctional(bw)
        summed_ts += -log10(compatibility(lut, f, p))
    end
    return summed_ts
end

"""
    function ECCDF(λ::Float64, lut::MinkowskiPValueLookup, N::Int, n::Int)

Calculates the ECCDF for a system size given by `lut` at average background `λ`,
using the LUT for fast p-value computation instead of building `MinkowskiDistribution`s.
"""
function ECCDF(λ::Float64, lut::MinkowskiPValueLookup, N::Int)
    L = lut.n
    ts = zeros(N)
    Threads.@threads for i in 1:N
        ts[i] = calc_ts(lut, CountsMap(L, λ), λ)
    end
    e_cdf = ecdf(ts)
    return ECCDF(λ, L, N, e_cdf)
end

"""
    function write_eccdf(path::AbstractString, x::ECCDF)

Stores the ECCDF as an HDF5 file at `path/eccdf_lambda=λ_L=L.h5`.
The sorted sample values of the underlying empirical CDF are written so that
the distribution can be reconstructed exactly on load.
"""
function write_eccdf(path::AbstractString, x::ECCDF)
    mkpath(path)
    h5open(joinpath(path, "eccdf_lambda=$(x.λ)_L=$(x.L).h5"), "w") do h5f
        write(h5f, "λ", x.λ)
        write(h5f, "L", x.L)
        write(h5f, "N", x.N)
        write(h5f, "sorted_values", x.F.sorted_values)
    end
end

"""
    function read_eccdf(fname::AbstractString)

Reads an ECCDF from an HDF5 file previously written by `write_eccdf`.
"""
function read_eccdf(fname::AbstractString)
    h5open(fname, "r") do h5f
        λ             = read(h5f, "λ")
        L             = read(h5f, "L")
        N             = read(h5f, "N")
        sorted_values = read(h5f, "sorted_values")
        ECCDF(λ, L, N, ecdf(sorted_values))
    end
end

"""
    function compatibility(eccdf::ECCDF, lut::MinkowskiPValueLookup, x::Union{CountsMap, Matrix{Int64}}, λ::Float64)

Returns the p-value for a given counts map by computing the test statistic via `lut`
and looking it up in `eccdf`.
"""
function compatibility(eccdf::ECCDF, lut::MinkowskiPValueLookup, x::Union{CountsMap, Matrix{Int64}}, λ::Float64)
    ts = calc_ts(lut, x, λ)
    return eccdf(ts)
end

"""
    function MinkowskiMap(x::CountsMap, b::Float64, lut::MinkowskiPValueLookup, eccdf::ECCDF)

Calculates a MinkowskiMap for a homogeneous background `b` using a precomputed
`MinkowskiPValueLookup` and `ECCDF`. The ECCDF must have been calibrated at the
same background value.
"""
function MinkowskiMap(x::CountsMap, b::Float64, lut::MinkowskiPValueLookup, eccdf::ECCDF)
    L = lut.n
    m, n = size(x.pixels)
    if isodd(L)
        l1 = floor(Int, L/2)
        l2 = floor(Int, L/2)
        αs    = zeros(n - 2l1, m - 2l2)
        signs = zeros(n - 2l1, m - 2l2)
    else
        l1 = floor(Int, (L-1)/2)
        l2 = ceil(Int, (L-1)/2)
        l  = l1 + l2
        αs    = zeros(n - l, m - l)
        signs = zeros(n - l, m - l)
    end
    Threads.@threads for j in l1+1:m-l2
        for i in l1+1:n-l2
            local_counts = x[i-l1:i+l2, j-l1:j+l2]
            pvalue = compatibility(eccdf, lut, local_counts, b)
            αs[i-l1, j-l1]    = pvalue
            signs[i-l1, j-l1] = mean(local_counts) > b ? 1.0 : -1.0
        end
    end
    return MinkowskiMap(αs .* signs)
end


"""
    function MinkowskiMap(x::CountsMap, b::Background, lut::MinkowskiPValueLookup, N::Int, n::Int)

Calculates a MinkowskiMap for an inhomogeneous background by computing the ECCDF
on-the-fly per unique λ value. Only one ECCDF is held in RAM at a time, making
this suitable when there are many unique background values. `N` is the number of
samples used to build each ECCDF, `n` is the number of discretization points stored.
"""
function MinkowskiMap(x::CountsMap, b::Background, lut::MinkowskiPValueLookup, N::Int)
    L = lut.n
    m, n_pix = size(x.pixels)
    if isodd(L)
        l1 = floor(Int, L/2)
        l2 = floor(Int, L/2)
        αs    = zeros(n_pix - 2l1, m - 2l2)
        signs = zeros(n_pix - 2l1, m - 2l2)
    else
        l1 = floor(Int, (L-1)/2)
        l2 = ceil(Int, (L-1)/2)
        l  = l1 + l2
        αs    = zeros(n_pix - l, m - l)
        signs = zeros(n_pix - l, m - l)
    end
    idx_dict = get_λ_idxs(b, L)
    Threads.@threads for λ in get_λs(b, L)
        λ == 0.0 && continue
        eccdf = ECCDF(λ, lut, N)
        for idx in idx_dict[λ]
            i, j = idx.I
            l1 < i <= m - l2 && l1 < j <= n_pix - l2 || continue
            local_counts     = x[i-l1:i+l2, j-l1:j+l2]
            local_background = b[i-l1:i+l2, j-l1:j+l2]
            correction!(local_counts, local_background, λ)
            pvalue = compatibility(eccdf, lut, local_counts, λ)
            αs[i-l1, j-l1]    = pvalue
            signs[i-l1, j-l1] = mean(local_counts) > λ ? 1.0 : -1.0
        end
    end
    return MinkowskiMap(αs .* signs)
end
