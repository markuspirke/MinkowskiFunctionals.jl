"""
    struct ECCDF

Datatype which stores the empirical cumulative distibution function
for summed p-values at different thresholds.
"""
struct ECCDF
    λ::Float64
    L::Int64
    N::Int64
    ts::Vector{Float64}
    pvalues::Vector{Float64}
end

"""
    function ECCDF(λ::Float64, L::Int64, ecdf::T, N) where {T <: ECDF}

Given an already calculated ecdf (from stats base) return the ECCDF.
"""
function ECCDF(λ::Float64, L::Int64, ecdf::T, N) where {T <: ECDF}
    xs = vcat(range(minimum(ecdf), maximum(ecdf), N))
    ys = 1 .- ecdf.(xs)

    ECCDF(λ, L, length(ecdf.sorted_values), xs, ys)
end

"""
    function ECCDF(ds_pvalues::Dict{Int64, Dict{MinkowskiFunctional, Float64}}, λ::Float64, L::Int64, N::Int64, n::Int64)

Calculates the empiriacl cumlative distribution function, given a dictory of a dictionary of functionals.
"""
function ECCDF(ds_pvalues::Dict{Int64, Dict{MinkowskiFunctional, Float64}}, λ::Float64, L::Int64, N::Int64, n::Int64)
    ts = zeros(N)
    Threads.@threads for i in 1:N
        ts[i] = calc_ts(ds_pvalues, CountsMap(L, λ))
    end
    e_cdf = ecdf(ts)
    xs = range(minimum(ts), maximum(ts), n)
    ys = 1.0 .- e_cdf.(xs)

    return ECCDF(λ, L, N, xs, ys)
end

"""
    function (eccdf::ECCDF)(x::Float64)

Finds the nearest pvalue within the ecdf. A self implemented
binary search algorithm is used, which converges quickly.
"""
function (eccdf::ECCDF)(x::Float64)
    idx = binary_search(eccdf.ts, x)
    return idx < length(eccdf.ts) ? eccdf.pvalues[idx] : 1/eccdf.N
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
    function write_eccdf(path::AbstractString, x::ECCDF; all_functionals=true)

Stores the ECCDF as a h5 file.
"""
function write_eccdf(path::AbstractString, x::ECCDF; all_functionals=true)
    if all_functionals
        h5open(joinpath(path, "eccdf_lambda=$(x.λ).h5"), "w") do h5f
            write(h5f, "λ", x.λ)
            write(h5f, "L", x.L)
            write(h5f, "N", x.N)
            write(h5f, "ts", x.ts)
            write(h5f, "pvalues", x.pvalues)
        end
    else
        h5open(joinpath(path, "eccdf_A_lambda=$(x.λ).h5"), "w") do h5f
            write(h5f, "λ", x.λ)
            write(h5f, "L", x.L)
            write(h5f, "N", x.N)
            write(h5f, "ts", x.ts)
            write(h5f, "pvalues", x.pvalues)
        end
    end
end
# is this needed?
function _check_eccdf_exists(path, λ)
    existing_pvalues = readdir(path)
    x = "eccdf_lambda=$(λ).ht"
    return x in existing_pvalues
end

"""
    function read_eccdf(fname::AbstractString)

Reads the ECCDF from a h5 file.
"""
function read_eccdf(fname::AbstractString)
    f = h5open(fname, "r")
    eccdf = ECCDF(read(f["λ"]), read(f["L"]), read(f["N"]), read(f["ts"]), read(f["pvalues"]))
    close(f)

    return eccdf
end

"""
    function ECCDF(λ::Float64, L::Int64, N, n)

Calculates the ecdf for a system size of L at an average background λ with the Area functional.
"""
function ECCDF(λ::Float64, L::Int64, N, n)
    ρmax = find_max_threshold(λ, L)
    d = Dict(ρ => AreaDistribution(L^2, λ, ρ) for ρ in 1:ρmax)
    ts = [calc_ts(d, CountsMap(L, λ)) for _ in 1:N]
    e_cdf = ecdf(ts)
    eccdf = ECCDF(λ, L, e_cdf, n)

    return eccdf
end

"""
    function compatibility(eccdf::ECCDF, d::Dict{Int64, T}, x::Union{CountsMap, Matrix{Int64}}) where {T<:AbstractMinkowskiDistribution}

Return the p-value for a given counts map or matrix.
"""
function compatibility(eccdf::ECCDF, d::Dict{Int64, T}, x::Union{CountsMap, Matrix{Int64}}) where {T<:AbstractMinkowskiDistribution}
    ts = calc_ts(d, x)

    return eccdf(ts)
end

"""
    function compatibility(eccdf::ECCDF, dd::DefaultDict{Int64, T, Int64}, x::Union{CountsMap, Matrix{Int64}}) where {T<:AbstractMinkowskiDistribution}

Return the p-value for a given counts map or matrix.
"""
function compatibility(eccdf::ECCDF, dd::DefaultDict{Int64, T, Int64}, x::Union{CountsMap, Matrix{Int64}}) where {T<:AbstractMinkowskiDistribution}
    ts = calc_ts(dd, x)

    return eccdf(ts)
end

"""
    function compatibility(eccdf::ECCDF, d::Dict{Int64, Dict{MinkowskiFunctional, Float64}}, x::Union{CountsMap, Matrix{Int64}})

Return the p-value for a given counts map or matrix.
"""
function compatibility(eccdf::ECCDF, d::Dict{Int64, Dict{MinkowskiFunctional, Float64}}, x::Union{CountsMap, Matrix{Int64}})
    ts = calc_ts(d, x)

    return eccdf(ts)
end


"""
    function compatibility(e_cdf::T, dd::DefaultDict{Int64, S, Int64}, x::Union{CountsMap, Matrix{Int64}}) where {T <: ECDF, S <: AbstractMinkowskiDistribution}

Return the p-value for a given counts map or matrix.
"""
function compatibility(e_cdf::T, dd::DefaultDict{Int64, S, Int64}, x::Union{CountsMap, Matrix{Int64}}) where {T <: ECDF, S <: AbstractMinkowskiDistribution}
    ts = calc_ts(dd, x)
    if 1.0 - e_cdf(ts) > 0.0
        return 1.0 - e_cdf(ts)
    else
        return 1/length(e_cdf.sorted_values)
    end
end

"""
    function calc_ts(d::Dict{Int64, Dict{MinkowskiFunctional, Float64}}, x::Union{CountsMap, Matrix{Int64}})

Calculates the test statistic for a given counts map or matrix.
"""
function calc_ts(d::Dict{Int64, Dict{MinkowskiFunctional, Float64}}, x::Union{CountsMap, Matrix{Int64}})
    ρs = get_thresholds(x)
    summed_ts = 0.0
    for ρ in ρs
        summed_ts += -log10(compatibility(d[ρ], ρ, x))
    end

    return summed_ts
end

"""
    function calc_ts(d::Dict{Int64, T}, x::Union{CountsMap, Matrix{Int64}}) where {T<:AbstractMinkowskiDistribution}

Calculates the test statistic for a given counts map or matrix.
"""
function calc_ts(d::Dict{Int64, T}, x::Union{CountsMap, Matrix{Int64}}) where {T<:AbstractMinkowskiDistribution}
    ρs = get_thresholds(x)
    summed_ts = 0.0
    for ρ in ρs
        summed_ts += -log10(compatibility(d[ρ], x))
    end

    return summed_ts
end

"""
    function calc_ts(dd::DefaultDict{Int64, T, Int64}, x::Union{CountsMap, Matrix{Int64}}) where {T<:AbstractMinkowskiDistribution}

Calculates the test statistic for a given counts map or matrix.
"""
function calc_ts(dd::DefaultDict{Int64, T, Int64}, x::Union{CountsMap, Matrix{Int64}}) where {T<:AbstractMinkowskiDistribution}
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

"""
    function MinkowskiMap(x::CountsMap, mink_ds::DefaultDict{Int64, S, Int64}, eccdf::T) where {S<:AbstractMinkowskiDistribution, T<:ECDF}

Calculates a minkowski map for a homogenous background, given a Dictionary of p-values and
a ECCDF.
"""
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
            pvalue = compatibility(eccdf, mink_ds, local_counts)
            αs[i-l, j-l] = pvalue
            signs[i-l, j-l] = mean(local_counts) > λ ? 1.0 : -1.0
        end
    end
    MinkowskiMap(αs .* signs)
end

"""
    function MinkowskiMap(x::CountsMap, mink_ds::Dict{Int64, Dict{MinkowskiFunctional, Float64}}, eccdf::ECCDF)

Calculates a minkowski map for a homogenous background, given a Dictionary of p-values and
a ECCDF.
"""
function MinkowskiMap(x::CountsMap, mink_ds::Dict{Int64, Dict{MinkowskiFunctional, Float64}}, eccdf::ECCDF)
    λ = eccdf.λ
    m, n = size(x)
    L = eccdf.L
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



"""
    function MinkowskiMap(x::CountsMap, mink_ds::Dict{Int64, Dict{MinkowskiFunctional, Float64}}, eccdf::ECCDF)

Calculates a minkowski map for a homogenous background based on the area functional, given a Dictionary of p-values and
a ECCDF.
"""
function MinkowskiMap(x::CountsMap, eccdf::ECCDF)
    ρs = get_thresholds(x)
    L = eccdf.L
    λ = eccdf.λ
    mink_ds = Dict(ρ => AreaDistribution(L^2, λ, ρ) for ρ in ρs)
    m, n = size(x)
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
