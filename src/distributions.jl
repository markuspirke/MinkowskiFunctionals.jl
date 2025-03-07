"""
    struct AreaDistribution

This defines a datatype which stores the area probability distribution
for a given system size and given λ and ρ.
"""
struct AreaDistribution
    n::Int
    λ::Float64
    ρ::Int
    p::Binomial{Float64}
end

function AreaDistribution(n, λ, ρ)
    d_poisson = Distributions.Poisson(λ)
    p = 1 - cdf(d_poisson, ρ-1)
    d_A = Binomial(n, p)

    AreaDistribution(n, λ, ρ, d_A)
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
struct MinkowskiDistribution
    n::Int
    p_black::Float64
    λ
    ρ::Int
    p::Accumulator
end

"""
    function MinkowskiDistribution(Ω::DensityOfStates, λ, ρ)

This generate the joint probability distribution out of the density of states.
"""
function MinkowskiDistribution(Ω::DensityOfStates, λ, ρ)
    p = 1 - cdf(Distributions.Poisson(λ), ρ-1)

    distribution = Accumulator{MinkowskiFunctional, Float64}()

    for (key, value) in Ω.data
        distribution[key] += value * p^key.A * (1 - p)^(Ω.n^2 - key.A)
    end

    MinkowskiDistribution(Ω.n, p, λ, ρ, distribution)#, c_distribution)
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
    p = pdf(d, f)
    ps = pdf(d)
    return sum(ps[ps .<= p])
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

# compatibility(d::MinkowskiDistribution, f::MinkowskiFunctional) = d.α[f]

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
