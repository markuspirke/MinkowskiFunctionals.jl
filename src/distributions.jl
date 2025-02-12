"""
    struct AreaDistribution

This defines a datatype which stores the area probability distribution
for a given system size and given λ and ρ.
"""
struct AreaDistribution
    n::Int
    λ
    ρ::Int
    p::Binomial{Float64}
    α::Dict{Int64, Float64}
end

function AreaDistribution(n, λ, ρ)
    d_poisson = Distributions.Poisson(λ)
    p = 1 - cdf(d_poisson, ρ-1)
    d_A = Binomial(n^2, p)
    α = get_αs_binomial(d_A)

    AreaDistribution(n, λ, ρ, d_A, α)
end

"""
    function get_αs_binomial(d::Binomial{Float64})

A small helper function to precalculate a lookup table for
hypothesis testing based on a Binomial distribution.
"""
function get_αs_binomial(d::Binomial{Float64})
    ps = pdf(d)
    αs = Dict(i-1 => sum(ps[ps[i] .>= ps]) for i in 1:length(ps))

    return αs
end

"""
    struct MinkowskiDistribution

This defines a datatype which stores the joint probability distribution
for a given system size and given λ and ρ.
"""
struct MinkowskiDistribution
    n::Int
    λ
    ρ::Int
    p::Accumulator
    α::Accumulator
end

function MinkowskiDistribution(n, λ, ρ, distribution::Accumulator)
    c_distribution =  deepcopy(distribution)
    ps = collect(values(distribution))
    for (key,value) in distribution # THIS IS NOT YET PERFECT, BECAUSE I DO NOT WANT TO USE IT
        p = distribution[key]
        mask = ps .<= p
        c_distribution[key] = sum(ps[mask])
    end

    MinkowskiDistribution(n, λ, ρ, distribution, c_distribution)
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

    c_distribution =  deepcopy(distribution)
    ps = collect(values(distribution))
    ps_perm = sortperm(ps) # THIS IS PROVISIONALLY
    ps = ps[ps_perm]
    c_ps = cumsum(ps)
    new_ps = c_ps[invperm(ps_perm)]
    for (i, key) in enumerate(keys(c_distribution))
        c_distribution[key] = new_ps[i]
    end
    # @showprogress for (key,p) in distribution
    #     # p = distribution[key]
    #     mask = ps .<= p
    #     c_distribution[key] = sum(ps[mask])
    # end

    MinkowskiDistribution(Ω.n, λ, ρ, distribution, c_distribution)
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

compatibility(d::MinkowskiDistribution, f::MinkowskiFunctional) = d.α[f]


function append!(h5f::HDF5.File, distribution::MinkowskiDistribution)
    n = distribution.n
    λ = distribution.λ
    ρ = distribution.ρ

    if n in parse.(Int, keys(h5f))
        # println("Macrostates for system size $(n) already exist.")
        if λ in parse.(Int, getindex.(split.(filter(x -> occursin("=", x), keys(h5f["$(n)"])), "="), 2))
            # println("Background already exists.")
            if ρ in parse.(Int, getindex.(split.(keys(h5f["$(n)/λ=$(λ)"]), "="), 2))
                # println("Treshold already exists. No changes are applied.")
            else
                ps = collect(values(distribution.p))
                αs = collect(values(distribution.α))
                xs = [(p, α) for (p, α) in zip(ps, αs)]
                write_dataset(h5f, "$(n)/λ=$(λ)/ρ=$(ρ)", xs)
                attributes(h5f["$(n)/λ=$(λ)/ρ=$(ρ)"])[Dates.format(now(), dateformat"yyyy-mm-dd")] = "v0.4.0"#Pkg.TOML.parse(read("Project.toml", String))["version"]
            end
        else
            ps = collect(values(distribution.p))
            αs = collect(values(distribution.α))
            xs = [(p, α) for (p, α) in zip(ps, αs)]
            write_dataset(h5f, "$(n)/λ=$(λ)/ρ=$(ρ)", xs)
            attributes(h5f["$(n)/λ=$(λ)/ρ=$(ρ)"])[Dates.format(now(), dateformat"yyyy-mm-dd")] = "v0.4.0"#Pkg.TOML.parse(read("Project.toml", String))["version"]
        end
    else
        # println("Macrostates do not yet exist. Generating them.")
        write_dataset(h5f, "$(n)/macrostates", collect(keys(distribution.α)))
        ps = collect(values(distribution.p))
        αs = collect(values(distribution.α))
        xs = [(p, α) for (p, α) in zip(ps, αs)]

        write_dataset(h5f, "$(n)/λ=$(λ)/ρ=$(ρ)", xs)
        attributes(h5f["$(n)/λ=$(λ)/ρ=$(ρ)"])[Dates.format(now(), dateformat"yyyy-mm-dd")] = "v0.4.0"#Pkg.TOML.parse(read("Project.toml", String))["version"]
    end
end
#     xs = @NamedTuple{A::Int64, P::Int64, χ::Int64, D::Float64, σ::Float64}[]
#     field_names = (:A, :P, :χ, :D, :α)
#     for k in collect(keys(distribution.P))
#         push!(xs, NamedTuple{field_names}((k.A, k.P, k.χ, distribution.P[k], distribution.σ[k])))
#     end
#     write_dataset(h5f, "$(n)/λ=$(λ)/ρ=$(ρ)", xs)
#     attributes(h5f["$(n)/λ=$(λ)/ρ=$(ρ)"])[Dates.format(now(), dateformat"yyyy-mm-dd")] = "v0.4.0"#Pkg.TOML.parse(read("Project.toml", String))["version"]
# end

function MinkowskiDistribution(fname::AbstractString, n::Int64, λ::Int64, ρ::Int64)
    h5open(fname, "r") do h5f
        d_counter = Accumulator{MinkowskiFunctional, Float64}()
        α_counter = Accumulator{MinkowskiFunctional, Float64}()

        functionals = read(h5f["$(n)/macrostates"])
        xs = read(h5f["$n/λ=$λ/ρ=$ρ"])
        for (f, x) in zip(functionals, xs)
            d_counter[MinkowskiFunctional(f.A, f.P, f.χ)] = x[1]
            α_counter[MinkowskiFunctional(f.A, f.P, f.χ)] = x[2]
        end

        MinkowskiDistribution(n, λ, ρ, d_counter, α_counter)
    end
end

function MinkowskiDistribution(h5f::HDF5.File, n::Int64, λ::Int64, ρ::Int64)
    d_counter = Accumulator{MinkowskiFunctional, Float64}()
    α_counter = Accumulator{MinkowskiFunctional, Float64}()

    functionals = read(h5f["$(n)/macrostates"])
    xs = read(h5f["$n/λ=$λ/ρ=$ρ"])
    for (f, x) in zip(functionals, xs)
        d_counter[MinkowskiFunctional(f.A, f.P, f.χ)] = x[1]
        α_counter[MinkowskiFunctional(f.A, f.P, f.χ)] = x[2]
    end

    MinkowskiDistribution(n, λ, ρ, d_counter, α_counter)
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
