"""
    struct MinkowskiDistribution

This defines a datatype which stores the joint probability distribution
for a given system size and given λ and ρ.
"""
struct MinkowskiDistribution
    n::Int
    λ::Int
    ρ::Int
    P::Accumulator
    σ::Accumulator

end

function MinkowskiDistribution(n, λ, ρ, distribution::Accumulator)
    c_distribution =  deepcopy(distribution)
    ps = collect(values(distribution))
    for (key,value) in distribution # THIS IS NOT YET PERFECT, BECAUSE I DO NOT WANT TO USE IT
        p = distribution[key]
        mask = ps .<= p
        c_distribution[key] = p2σ(sum(ps[mask]))
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
    new_ps = p2σ.(c_ps[invperm(ps_perm)])
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
    return d.P[f]
end

function Distributions.pdf(d::MinkowskiDistribution)
    return collect(values(d.P))
end

compatibility(d::MinkowskiDistribution, f::MinkowskiFunctional) = d.σ[f]


function append!(h5f::HDF5.File, distribution::MinkowskiDistribution)
    n = distribution.n
    λ = distribution.λ
    ρ = distribution.ρ

    if n in parse.(Int, keys(h5f))
        println("Macrostates for system size $(n) already exist.")
        if λ in parse.(Int, getindex.(split.(filter(x -> occursin("=", x), keys(h5f["$(n)"])), "="), 2))
            println("Background already exists.")
            if ρ in parse.(Int, getindex.(split.(keys(h5f["$(n)/λ=$(λ)"]), "="), 2))
                println("Treshold already exists. No changes are applied.")
            else
                ps = collect(values(distribution.P))
                σs = collect(values(distribution.σ))
                xs = [(p, σ) for (p, σ) in zip(ps, σs)]
                write_dataset(h5f, "$(n)/λ=$(λ)/ρ=$(ρ)", xs)
                attributes(h5f["$(n)/λ=$(λ)/ρ=$(ρ)"])[Dates.format(now(), dateformat"yyyy-mm-dd")] = "v0.4.0"#Pkg.TOML.parse(read("Project.toml", String))["version"]
            end
        else
            ps = collect(values(distribution.P))
            σs = collect(values(distribution.σ))
            xs = [(p, σ) for (p, σ) in zip(ps, σs)]
            write_dataset(h5f, "$(n)/λ=$(λ)/ρ=$(ρ)", xs)
            attributes(h5f["$(n)/λ=$(λ)/ρ=$(ρ)"])[Dates.format(now(), dateformat"yyyy-mm-dd")] = "v0.4.0"#Pkg.TOML.parse(read("Project.toml", String))["version"]
        end
    else
        println("Macrostates do not yet exist. Generating them.")
        write_dataset(h5f, "$(n)/macrostates", collect(keys(distribution.σ)))
        ps = collect(values(distribution.P))
        σs = collect(values(distribution.σ))
        xs = [(p, σ) for (p, σ) in zip(ps, σs)]

        write_dataset(h5f, "$(n)/λ=$(λ)/ρ=$(ρ)", xs)
        attributes(h5f["$(n)/λ=$(λ)/ρ=$(ρ)"])[Dates.format(now(), dateformat"yyyy-mm-dd")] = "v0.4.0"#Pkg.TOML.parse(read("Project.toml", String))["version"]
    end
end
#     xs = @NamedTuple{A::Int64, P::Int64, χ::Int64, D::Float64, σ::Float64}[]
#     field_names = (:A, :P, :χ, :D, :σ)
#     for k in collect(keys(distribution.P))
#         push!(xs, NamedTuple{field_names}((k.A, k.P, k.χ, distribution.P[k], distribution.σ[k])))
#     end
#     write_dataset(h5f, "$(n)/λ=$(λ)/ρ=$(ρ)", xs)
#     attributes(h5f["$(n)/λ=$(λ)/ρ=$(ρ)"])[Dates.format(now(), dateformat"yyyy-mm-dd")] = "v0.4.0"#Pkg.TOML.parse(read("Project.toml", String))["version"]
# end

function MinkowskiDistribution(fname::AbstractString, n::Int64, λ::Int64, ρ::Int64)
    h5open(fname, "r") do h5f
        d_counter = Accumulator{MinkowskiFunctional, Float64}()
        σ_counter = Accumulator{MinkowskiFunctional, Float64}()

        functionals = read(h5f["$(n)/macrostates"])
        xs = read(h5f["$n/λ=$λ/ρ=$ρ"])
        for (f, x) in zip(functionals, xs)
            d_counter[MinkowskiFunctional(f.A, f.P, f.χ)] = x[1]
            σ_counter[MinkowskiFunctional(f.A, f.P, f.χ)] = x[2]
        end

        MinkowskiDistribution(n, λ, ρ, d_counter, σ_counter)
    end
end

function marginalize(P::MinkowskiDistribution, field)
    marginalized_distribution = Accumulator{Int64, Float64}()

    for (x, p) in P.P
        marginalized_distribution[getfield(x, field)] += p
    end


    DiscreteNonParametric(collect(keys(marginalized_distribution)), collect(values(marginalized_distribution)))
end
