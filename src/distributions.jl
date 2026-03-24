abstract type AbstractMinkowskiDistribution end
"""
    struct AreaDistribution <: AbstractMinkowskiDistribution

This defines a datatype which stores the area probability distribution
for a given system size and given λ and ρ.
"""
struct AreaDistribution <: AbstractMinkowskiDistribution
    n::Int
    λ::Union{Float64, Vector{Float64}}
    ρ::Int
    p::Union{Binomial{Float64}, PoissonBinomial{Float64, Vector{Float64}}}
    # p::Binomial{Float64}
    pvalues::Union{Accumulator, Missing}
end

function Statistics.mean(d::AreaDistribution)
    return mean(d.p)
end

"""
    function window_size(d::AreaDistribution)

Returns the system size of the given distribution.
"""
function window_size(d::AreaDistribution)
    return Int(sqrt(d.n))
end


"""
    function AreaDistribution(n, λ, ρ; pvalues=true)

This calculates the Area distribution for a given system size, background and threshold.
By default is also precalculates the pvalues.
"""
function AreaDistribution(n::Int, λ::Float64, ρ::Int; pvalues=true)
    p, _ = gamma_inc(ρ, λ)
    d_A = Binomial(n, p)

    if pvalues
        d_pvalues = get_pvalues(d_A)
        return AreaDistribution(n, λ, ρ, d_A, d_pvalues)
    else
        return AreaDistribution(n, λ, ρ, d_A, missing)
    end
end

"""
    function AreaDistribution(λs::Matrix{Float64}, ρ; pvalues=true)

This calculates the Area distribution for a given, possibly nonhomgenous background and threshold.
By default is also precalculates the pvalues.
"""
function AreaDistribution(λs::Matrix{Float64}, ρ; pvalues=true)
    ps = [gamma_inc(ρ, λ)[1] for λ in λs[:] if λ > 0.0]
    d_A = PoissonBinomial(ps)

    if pvalues
        d_pvalues = get_pvalues(d_A)
        return AreaDistribution(length(ps), λs[:], ρ, d_A, d_pvalues)
    else
        return AreaDistribution(length(ps), λs[:], ρ, d_A, missing)
    end
end

function AreaDistribution(λs::Background, ρ; pvalues=true)
    AreaDistribution(λs.pixels, ρ; pvalues=pvalues)
end

function get_pvalues(d::Binomial{Float64})
    ks, vs = support(d), pdf(d)
    idxs = sortperm(vs)
    ks, vs = ks[idxs], vs[idxs]
    ps = zeros(length(ks))
    modified_cumsum!(ps, vs)
    d_pvalues = Dict(k => p for (k, p) in zip(ks, ps))

    return Accumulator(d_pvalues)

end

function get_pvalues(d::PoissonBinomial{Float64, Vector{Float64}})
    ks, vs = support(d), pdf(d)
    idxs = sortperm(vs)
    ks, vs = ks[idxs], vs[idxs]
    ps = zeros(length(ks))
    modified_cumsum!(ps, vs)
    d_pvalues = Dict(k => p for (k, p) in zip(ks, ps))

    return Accumulator(d_pvalues)

end

function Base.show(io::IO, P::AreaDistribution)
    print(io, "Area distribution for n=$(P.n), λ=$(P.λ) and ρ=$(P.ρ).")
end

"""
    function Distributions.pdf(d::AreaDistribution, f::MinkowskiFunctional)

Returns the probability for a given functional with respect to the Area distribution.
"""
function Distributions.pdf(d::AreaDistribution, f::MinkowskiFunctional)
    return pdf(d.p, f.A)
end

"""
    function Distributions.pdf(d::AreaDistribution, f::T) where {T<:Integer}

Returns the probaility for a given area value.
"""
function Distributions.pdf(d::AreaDistribution, f::T) where {T<:Integer}
    return pdf(d.p, f)
end

"""
    function Distributions.pdf(d::AreaDistribution)

Returns the probability density function as a 1D vector.
"""
function Distributions.pdf(d::AreaDistribution)
    return pdf(d.p)
end

"""
    function compatibility(d::AreaDistribution, f::MinkowskiFunctional)

Returns the p-values for a given functional.
"""
function compatibility(d::AreaDistribution, f::MinkowskiFunctional)
    return compatibility(d, f.A)
end

"""
    function compatibility(d::AreaDistribution, f::T) where {T<:Integer}

Return the p-value for a given area A.
"""
function compatibility(d::AreaDistribution, f::T) where {T<:Integer}
    if d.pvalues |> ismissing
        p = pdf(d, f)
        ps = pdf(d)
        return sum(ps[ps .<= p])
    else
        return d.pvalues[f]
    end
end

function compatibility(d::AreaDistribution, x::CountsMap)
    bw_map = BWMap(x, d.ρ)
    A = sum(bw_map.pixels)
    return compatibility(d, A)
end

function compatibility(d::AreaDistribution, x::Matrix{Int64})
    bw_map = BWMap(x, d.ρ)
    A = sum(bw_map.pixels)
    return compatibility(d, A)
end

"""
    struct MinkowskiDistribution

This defines a datatype which stores the joint probability distribution
for a given system size and given λ and ρ.
"""
struct MinkowskiDistribution <: AbstractMinkowskiDistribution
    n::Int
    p_black::Float64
    λ
    ρ::Int
    p::Accumulator
    pvalues::Union{Accumulator, Missing}
end


"""
    function window_size(d::MinkowskiDistribution)

Returns the system size of the given distribution.
"""
function window_size(d::MinkowskiDistribution)
    return d.n
end

"""
    function MinkowskiDistribution(Ω::DensityOfStates, λ, ρ)

This generate the joint probability distribution out of the density of states.
"""
function MinkowskiDistribution(Ω::DensityOfStates, λ, ρ; pvalues=true)
    # p = 1 - cdf(Distributions.Poisson(λ), ρ-1)
    p, _ = gamma_inc(ρ, λ)
    distribution = Accumulator{MinkowskiFunctional, Float64}()

    for (key, value) in Ω.data
        distribution[key] += value * p^key.A * (1 - p)^(Ω.n^2 - key.A)
    end
    mink_distribution = MinkowskiDistribution(Ω.n, p, λ, ρ, distribution, missing)

    if pvalues
        d_pvalues = get_pvalues(distribution)
        return mink_distribution = MinkowskiDistribution(Ω.n, p, λ, ρ, distribution, d_pvalues)
    else
        return mink_distribution
    end
end

function get_pvalues(d::Accumulator{MinkowskiFunctional, Float64})
    ks, vs = collect(d |> keys), collect(d |> values)
    idxs = sortperm(vs)
    ks, vs = ks[idxs], vs[idxs]
    ps = zeros(length(ks))
    modified_cumsum!(ps, vs)
    d_pvalues = Dict(k => p for (k, p) in zip(ks, ps))

    return Accumulator(d_pvalues)

end

function modified_cumsum!(A, B)
    x = 0.0
    i = 1
    while i <= length(B)
        b = B[i]
        j = i
        # Find how many times b repeats
        while j <= length(B) && B[j] == b
            j += 1
        end
        count = j - i
        total = b * count
        x += total
        for k in i:j-1
            A[k] = x
        end
        i = j
    end
end

function Base.show(io::IO, P::MinkowskiDistribution)
    print(io, "Minkowski distribution for n=$(P.n), λ=$(P.λ) and ρ=$(P.ρ).")
end
"""
    function Distributions.pdf(d::MinkowskiDistribution, f::MinkowskiFunctional)

Given a Minkowski distribution, this evaluates the probabilty of the given Minkowksi functional.
"""
function Distributions.pdf(d::MinkowskiDistribution, f::MinkowskiFunctional)
    return d.p[f]
end

"""
    function Distributions.pdf(d::MinkowskiDistribution, f::NamedTuple)

Given a Minkowski distribution, this evaluates the probabilty of the given Minkowksi functional.
"""
function Distributions.pdf(d::MinkowskiDistribution, f::NamedTuple)
    return d.p[f]
end

"""
    function Distributions.pdf(d::MinkowskiDistribution)

Given a Minkowski distribution, this return the PDF as a single one dimensional vector.
"""
function Distributions.pdf(d::MinkowskiDistribution)
    return collect(values(d.p))
end

"""
    function compatibility(d::MinkowskiDistribution, f::MinkowskiFunctional)

This returns the p-value for a given Minkowski functional.
"""
function compatibility(d::MinkowskiDistribution, f::MinkowskiFunctional)
    if d.pvalues |> ismissing
        p = pdf(d, f)
        ps = pdf(d)
        return sum(ps[ps .<= p])
    else
        return d.pvalues[f]
    end
end
"""
    function compatibility(d::MinkowskiDistribution, x::CountsMap)

This return the p-value of a given counts map at the threshold of the Minkowski distribution.
"""
function compatibility(d::MinkowskiDistribution, x::CountsMap)
    compatibility(d, x.pixels)
end

"""
    function compatibility(d::MinkowskiDistribution, x::CountsMap)

This return the p-value of a given counts map at the threshold of the Minkowski distribution.
"""
function compatibility(d::MinkowskiDistribution, x::Matrix{Int64})
    bw_map = BWMap(x, d.ρ)
    f = MinkowskiFunctional(bw_map)
    return compatibility(d, f)
end

"""
    function marginalize(P::MinkowskiDistribution, fields)

This takes the 3 dimensional Minkowski Distribution and marginalizes along all given fields.
It returns again a Minkowski distribution.
"""
function marginalize(P::MinkowskiDistribution, fields)
    old_fields = (:A, :P, :χ)
    new_fields = Symbol[]
    for field in old_fields
        if !(field in fields)
            push!(new_fields, field)
        end
    end
    new_fields = Tuple(new_fields)
    counter_P = Accumulator{NamedTuple{new_fields}, Float64}()

    for k in keys(P.p)
        k_reduced = reduce_functional(k, new_fields)
        counter_P[k_reduced] += P.p[k]
    end

    counter_α =  deepcopy(counter_P)
    ps = collect(values(counter_P))
    ps_perm = sortperm(ps) # THIS IS PROVISIONALLY
    ps = ps[ps_perm]
    c_ps = cumsum(ps)
    new_ps = c_ps[invperm(ps_perm)]
    for (i, key) in enumerate(keys(counter_α))
        counter_α[key] = new_ps[i]
    end


    return MinkowskiDistribution(P.n, P.p_black, P.λ, P.ρ, counter_P, counter_α)
end

function reduce_functional(functional::MinkowskiFunctional, fields::T) where T<:Tuple
    return NamedTuple{fields}(Tuple([getfield(functional, f) for f in fields]))
end


"""
    function MinkowskiDistribution(h5f::HDF5.File, n::Int64, λ::Float64, ρ::Int64)

Loads a Minkwoski distribution from a given h5 file.
"""
function MinkowskiDistribution(h5f::HDF5.File, n::Int64, λ::Float64, ρ::Int64)
    d_counter = Accumulator{MinkowskiFunctional, Float64}()
    pvalues_counter = Accumulator{MinkowskiFunctional, Float64}()
    p_black = 1 - cdf(Distributions.Poisson(λ), ρ-1)

    functionals = read(h5f["$(n)/macrostates"])
    ps = read(h5f["$n/λ=$λ/ρ=$ρ"])
    for (f, x) in zip(functionals, ps)
        d_counter[MinkowskiFunctional(f.A, f.P, f.χ)] = x[1]
        pvalues_counter[MinkowskiFunctional(f.A, f.P, f.χ)] = x[2]
    end
    if pvalues_counter[MinkowskiFunctional(0, 0, 0)] == -1
        MinkowskiDistribution(n, p_black, λ, ρ, d_counter, missing)
    else
        MinkowskiDistribution(n, p_black, λ, ρ, d_counter, pvalues_counter)
    end
end

"""
    function append!(h5f::HDF5.File, distribution::MinkowskiDistribution)

Appends a Minkowski distribution to a given h5 file.
"""
function append!(h5f::HDF5.File, distribution::MinkowskiDistribution)
    n = distribution.n
    λ = distribution.λ
    ρ = distribution.ρ
    if n in parse.(Int, keys(h5f)) # IF FILE HAS ENTRY FOR SYSTEM SIZE N
        println("Macrostates for system size $(n) already exist.")
        if λ in parse.(Float64, getindex.(split.(filter(x -> occursin("=", x), keys(h5f["$(n)"])), "="), 2)) # CHECK ENTRIES WITH SYSTEM SIZE n
            println("Background already exists.")
            if ρ in parse.(Int, getindex.(split.(keys(h5f["$(n)/λ=$(λ)"]), "="), 2)) # CHECK IF THRESHOLD
                println("Threshold already exists. No changes are applied.")
            else
                ps = collect(values(distribution.p))
                if ismissing(distribution.pvalue)
                    pvalues = [-1 for _ in 1:length(ps)] # IF NO PVALUES RETURN -1
                else
                    pvalues = collect(values(distribution.pvalues))
                end
                xs = [(p, α) for (p, α) in zip(ps, pvalues)]
                write_dataset(h5f, "$(n)/λ=$(λ)/ρ=$(ρ)", xs)
                # write_dataset(h5f, "$(n)/λ=$(λ)/ρ=$(ρ)", ps)
                # attributes(h5f["$(n)/λ=$(λ)/ρ=$(ρ)"])[Dates.format(now(), dateformat"yyyy-mm-dd")] = "v0.4.0"#Pkg.TOML.parse(read("Project.toml", String))["version"]
            end
        else
            ps = collect(values(distribution.p))
            if ismissing(distribution.pvalue)
                pvalues = [-1 for _ in 1:length(ps)]
            else
                pvalues = collect(values(distribution.pvalues))
            end
            xs = [(p, α) for (p, α) in zip(ps, pvalues)]
            write_dataset(h5f, "$(n)/λ=$(λ)/ρ=$(ρ)", xs)
            # ps = collect(values(distribution.p))
            # pvalues = collect(values(distribution.pvalue))
            # xs = [(p, α) for (p, α) in zip(ps, pvalues)]
            # write_dataset(h5f, "$(n)/λ=$(λ)/ρ=$(ρ)", xs)
            # attributes(h5f["$(n)/λ=$(λ)/ρ=$(ρ)"])[Dates.format(now(), dateformat"yyyy-mm-dd")] = "v0.4.0"#Pkg.TOML.parse(read("Project.toml", String))["version"]
        end
    else
        println("Macrostates do not yet exist. Generating them.")

        write_dataset(h5f, "$(n)/macrostates", collect(keys(distribution.p)))

        ps = collect(values(distribution.p))
        if ismissing(distribution.pvalues)
            pvalues = [-1 for _ in 1:length(ps)]
        else
            pvalues = collect(values(distribution.pvalues))
        end
        xs = [(p, α) for (p, α) in zip(ps, pvalues)]
        write_dataset(h5f, "$(n)/λ=$(λ)/ρ=$(ρ)", xs)
        # ps = collect(values(distribution.p))
        # pvalues = collect(values(distribution.pvalue))
        # xs = [(p, α) for (p, α) in zip(ps, pvalues)]
        # write_dataset(h5f, "$(n)/λ=$(λ)/ρ=$(ρ)", xs)
        # attributes(h5f["$(n)/λ=$(λ)/ρ=$(ρ)"])[Dates.format(now(), dateformat"yyyy-mm-dd")] = "v0.4.0"#Pkg.TOML.parse(read("Project.toml", String))["version"]
    end
end

struct MinkowskiFunctionalX
    A::Int16
    P::Int16
    χ::Int16
    p::Float64
end


"""

This stores a Minkowski distribution as a binary file at the given path with
a filename that is like lambda=λ_rho=ρ.dat.
"""
function write_pvalues(path::AbstractString, d::MinkowskiDistribution)
    xs = d.pvalues
    ys = [MinkowskiFunctionalX(k.A, k.P, k.χ, v) for (k, v) in xs]
    write(joinpath(path, "lambda=$(d.λ)_rho=$(d.ρ).dat"), ys)
end

"""

This stores all needed Minkowski distribution to sample the eccdf,
as a binary file at the given path with a filename that is like lambda=λ_rho=ρ.dat.
"""
function write_necessary_pvalues(path::AbstractString, b::Background, Ω::DensityOfStates)
    L = Ω.n
    λs = get_λs(b, L)
    ρs = [find_max_threshold(λ, L) for λ in λs]
    for (λ, ρmax) in zip(λs, ρs)
        for ρ in 1:ρmax
            _check_exists(path, λ, ρ) && continue
            write_pvalues(path, MinkowskiDistribution(Ω, λ, ρ))
        end
    end
end

function write_necessary_pvalues(path::AbstractString, λ::Float64, Ω::DensityOfStates)
    L = Ω.n
    ρmax = find_max_threshold(λ, L)
    for ρ in 1:ρmax
        _check_exists(path, λ, ρ) && continue
        write_pvalues(path, MinkowskiDistribution(Ω, λ, ρ))
    end
end
"""

This stores all needed Minkowski distribution to calculate a sky map,
as a binary file at the given path with a filename that is like lambda=λ_rho=ρ.dat.
"""
function write_necessary_pvalues(path::AbstractString, b::Background, x::CountsMap, Ω::DensityOfStates)
    L = Ω.n
    m, n = size(x)
    d_λ_idxs = MinkowskiFunctionals.get_λ_idxs(b, L)
    MinkowskiFunctionals.remove_boundary_idxs!(d_λ_idxs, m, n, L)
    d_ρ_λ = MinkowskiFunctionals.get_ρ_λ(x, d_λ_idxs, L)
    for (λ, ρs) in d_ρ_λ
        for ρ in ρs
            # @show λ, ρ, path
            _check_exists(path, λ, ρ) && continue
            write_pvalues(path, MinkowskiDistribution(Ω, λ, ρ))
        end
    end
end

"""
    function _check_exists(path, λ, ρ)

Helper function which checks whether p values are already calculated and stored at the
given path, λ and threshold.
"""
function _check_exists(path, λ, ρ)
    existing_pvalues = readdir(path)
    x = "lambda=$(λ)_rho=$(ρ).dat"
    return x in existing_pvalues
end

"""
    function read_pvalues(fname::AbstractString)

This reads stored p-values and return a dictory, which can be index by functionals.
"""
function read_pvalues(fname::AbstractString)
    xs = reinterpret(MinkowskiFunctionalX, read(fname))
    return Dict(MinkowskiFunctional(x.A, x.P, x.χ) => x.p for x in xs)
end

"""
    function read_pvalues(path::AbstractString, λ)

This reads all p-values for a given λ withing the same directory.
"""
function read_pvalues(path::AbstractString, λ)
    fnames = readdir(path, join=true)
    filter!(x -> occursin("lambda=$(λ)_rho", x), fnames)
    Dict(parse(Int64, match(r"rho=(\d+)\.", fname).captures[1]) => read_pvalues(fname) for fname in fnames)
end

"""
    function compatibility(d::Dict{MinkowskiFunctional, Float64}, ρ::Int64, x::Matrix{Int64})

This calculates the p-value for a given counts map (as a matrix).
"""
function compatibility(d::Dict{MinkowskiFunctional, Float64}, ρ::Int64, x::Matrix{Int64})
    bw_map = BWMap(x, ρ)
    f = MinkowskiFunctional(bw_map)
    return d[f]
end

"""
    function compatibility(d::Dict{MinkowskiFunctional, Float64}, ρ::Int64, x::CountsMap)

This calculates the p-value for a given counts map at the given threshold.
"""
function compatibility(d::Dict{MinkowskiFunctional, Float64}, ρ::Int64, x::CountsMap)
    compatibility(d, ρ, x.pixels)
end


"""
    struct MinkowskiPValueLookup

Precomputed lookup structure built once from a `DensityOfStates` for
allocation-free p-value queries at any blackening probability `p`.

Instead of building a full `MinkowskiDistribution` (an `Accumulator` dict) for
every `(λ, ρ)` pair, this struct stores each macrostate's Ω count in a sorted
array grouped by area value A, with precomputed cumulative sums.

At query time, the p-value of an observed `MinkowskiFunctional` `f` at blackening
probability `p = gamma_inc(ρ, λ)[1]` is computed in O(n² · log m) time —
where m is the max group size — using one `searchsortedlast` per distinct area
value and zero allocations.

Construction is O(N log N) in the number of macrostates and is done once per
`DensityOfStates`.
"""
struct MinkowskiPValueLookup
    n::Int
    n_sq::Int
    A_Ωs::Vector{Vector{Float64}}    # Ω counts, sorted ascending within each A-group
    A_cumΩ::Vector{Vector{Float64}}  # cumulative Ω sums within each A-group
    functional_Ω::Dict{MinkowskiFunctional, Float64}  # Ω lookup for observed f
end

"""
    function MinkowskiPValueLookup(Ω::DensityOfStates)

Builds a `MinkowskiPValueLookup` from a `DensityOfStates` by grouping all
macrostates by their area A and sorting within each group by Ω ascending.
"""
function MinkowskiPValueLookup(Ω::DensityOfStates)
    n    = Ω.n
    n_sq = n^2

    A_raw        = [Float64[] for _ in 0:n_sq]
    functional_Ω = Dict{MinkowskiFunctional, Float64}()

    for (f, count) in Ω.data
        A   = Int(f.A)
        Ω_f = Float64(count)
        push!(A_raw[A + 1], Ω_f)
        functional_Ω[f] = Ω_f
    end

    A_Ωs  = Vector{Vector{Float64}}(undef, n_sq + 1)
    A_cumΩ = Vector{Vector{Float64}}(undef, n_sq + 1)

    for A in 0:n_sq
        sorted       = sort(A_raw[A + 1])
        A_Ωs[A + 1]  = sorted
        A_cumΩ[A + 1] = cumsum(sorted)
    end

    MinkowskiPValueLookup(n, n_sq, A_Ωs, A_cumΩ, functional_Ω)
end

function Base.show(io::IO, lut::MinkowskiPValueLookup)
    print(io, "MinkowskiPValueLookup for a system size of $(lut.n) x $(lut.n).")
end

function window_size(lut::MinkowskiPValueLookup)
    return lut.n
end

"""
    function compatibility(lut::MinkowskiPValueLookup, f::MinkowskiFunctional, p::Float64)

Returns the p-value for an observed `MinkowskiFunctional` `f` at blackening
probability `p = gamma_inc(ρ, λ)[1]`.

The p-value is the sum of probabilities of all macrostates whose probability
does not exceed that of `f`. For macrostates in area group A, the inclusion
condition reduces to

    Ω(A, P, χ) ≤ Ω(f) · (p / (1 − p))^(f.A − A)

which is evaluated allocation-free via `searchsortedlast` in each sorted
per-area array.
"""
function compatibility(lut::MinkowskiPValueLookup, f::MinkowskiFunctional, p::Float64)
    Ω_obs = lut.functional_Ω[f]
    A_obs = Int(f.A)
    n_sq  = lut.n_sq

    p == 0.0 && return A_obs == 0    ? 1.0 : 0.0
    p == 1.0 && return A_obs == n_sq ? 1.0 : 0.0

    log_p   = log(p)
    log_1mp = log1p(-p)
    r       = p / (1 - p)

    pval = 0.0

    # Sweep downward: A = A_obs, A_obs-1, ..., 0.
    # threshold(A) = Ω_obs · r^(A_obs − A).  Starting at Ω_obs and multiplying by r
    # each step avoids any exp(log(Ω_obs)) round-trip, which would otherwise cause
    # searchsortedlast to miss the observed state in its own A-group.
    threshold = Ω_obs
    @inbounds for A in A_obs:-1:0
        group_Ωs = lut.A_Ωs[A + 1]
        if !isempty(group_Ωs)
            k = searchsortedlast(group_Ωs, threshold)
            if k > 0
                pval += lut.A_cumΩ[A + 1][k] * exp(A * log_p + (n_sq - A) * log_1mp)
            end
        end
        threshold *= r   # threshold for A − 1
    end

    # Sweep upward: A = A_obs+1, ..., n_sq.
    # threshold(A) = Ω_obs · r^(A_obs − A) = Ω_obs / r^(A − A_obs).
    threshold = Ω_obs / r
    @inbounds for A in A_obs+1:n_sq
        group_Ωs = lut.A_Ωs[A + 1]
        if !isempty(group_Ωs)
            k = searchsortedlast(group_Ωs, threshold)
            if k > 0
                pval += lut.A_cumΩ[A + 1][k] * exp(A * log_p + (n_sq - A) * log_1mp)
            end
        end
        threshold /= r   # threshold for A + 1
    end

    return pval
end

"""
    function sample_minkowski_distribution(b::Background, N::Int64)

This is preliminary. Was used to sample Functionals for a given non-uniform background.
"""
function sample_minkowski_distribution(b::Background, N::Int64)
    d_counter = Dict(ρ => Accumulator{MinkowskiFunctional, Int64}() for ρ in 1:100)
    for _ in 1:N
        counts_map = CountsMap(b)
        ρs = get_thresholds(counts_map)
        for ρ in ρs
            bw_map = BWMap(counts_map, ρ)
            f = MinkowskiFunctional(bw_map)
            d_counter[ρ][f] += 1
        end
    end
    d_counter
end
