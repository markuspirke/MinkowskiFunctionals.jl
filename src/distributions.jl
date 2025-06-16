abstract type AbstractMinkowskiDistribution end
"""
    struct AreaDistribution

This defines a datatype which stores the area probability distribution
for a given system size and given λ and ρ.
"""
struct AreaDistribution <: AbstractMinkowskiDistribution
    n::Int
    λ::Float64
    ρ::Int
    p::Binomial{Float64}
    pvalues::Union{Accumulator, Missing}
end

function window_size(d::AreaDistribution)
    return Int(sqrt(d.n))
end


function AreaDistribution(n, λ, ρ; pvalues=true)
    p, _ = gamma_inc(ρ, λ)
    d_A = Binomial(n, p)

    if pvalues
        d_pvalues = get_pvalues(d_A)
        return AreaDistribution(n, λ, ρ, d_A, d_pvalues)
    else
        return AreaDistribution(n, λ, ρ, d_A, missing)
    end
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

function Base.show(io::IO, P::AreaDistribution)
    print(io, "Area distribution for n=$(P.n), λ=$(P.λ) and ρ=$(P.ρ).")
end

function Distributions.pdf(d::AreaDistribution, f::MinkowskiFunctional)
    return pdf(d.p, f.A)
end

function Distributions.pdf(d::AreaDistribution, f::Int64)
    return pdf(d.p, f)
end

function Distributions.pdf(d::AreaDistribution)
    return pdf(d.p)
end

function compatibility(d::AreaDistribution, f::MinkowskiFunctional)
    p = pdf(d, f)
    ps = pdf(d)
    return sum(ps[ps .<= p])
end

function compatibility(d::AreaDistribution, f::Int64)
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
    f = MinkowskiFunctional(bw_map)
    return compatibility(d, f)
end

function compatibility(d::AreaDistribution, x::Matrix{Int64})
    bw_map = BWMap(x, d.ρ)
    f = MinkowskiFunctional(bw_map)
    return compatibility(d, f)
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

function Distributions.pdf(d::MinkowskiDistribution, f::MinkowskiFunctional)
    return d.p[f]
end

function Distributions.pdf(d::MinkowskiDistribution)
    return collect(values(d.p))
end


function compatibility(d::MinkowskiDistribution, f::MinkowskiFunctional)
    if d.pvalues |> ismissing
        p = pdf(d, f)
        ps = pdf(d)
        return sum(ps[ps .<= p])
    else
        return d.pvalues[f]
    end
end

function compatibility(d::MinkowskiDistribution, x::CountsMap)
    bw_map = BWMap(x, d.ρ)
    f = MinkowskiFunctional(bw_map)
    return compatibility(d, f)
end

function compatibility(d::MinkowskiDistribution, x::Matrix{Int64})
    bw_map = BWMap(x, d.ρ)
    f = MinkowskiFunctional(bw_map)
    return compatibility(d, f)
end


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


    return MinkowskiDistribution(P.n, P.λ, P.ρ, counter_P, counter_α)
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

function write_pvalues(path::AbstractString, d::MinkowskiDistribution)
    xs = d.pvalues
    ys = [MinkowskiFunctionalX(k.A, k.P, k.χ, v) for (k, v) in xs]
    write(joinpath(path, "lambda=$(d.λ)_rho=$(d.ρ).dat"), ys)
end

function read_pvalues(fname::AbstractString)
    xs = reinterpret(MinkowskiFunctionalX, read(fname))
    return Dict(MinkowskiFunctional(x.A, x.P, x.χ) => x.p for x in xs)
end

function compatibility(d::Dict{MinkowskiFunctional, Float64}, ρ::Int64, x::Matrix{Int64})
    bw_map = BWMap(x, ρ)
    f = MinkowskiFunctional(bw_map)
    return d[f]
end

function compatibility(d::Dict{MinkowskiFunctional, Float64}, ρ::Int64, x::CountsMap)
    compatibility(d, ρ, x.pixels)
end
