struct DensityOfStates
    n
    data::Accumulator{MinkowskiFunctional, Integer}
end

"""
    function DensityOfStates(n)

This calculates the Density of States for a given system size n.
"""
function DensityOfStates(n)
    counts = Accumulator{MinkowskiFunctional, Int64}()
    n_microstates = 2^(n^2)

    state = BitMatrix(undef, n, n)
    # Iterate over all 2^(n^2) binary configurations
    for i in 0:n_microstates-1
        # Fill the matrix `state` using bitwise operations
        x = i
        for k in 1:n^2
            state[k] = x & 1  # Set each element based on the binary representation
            x >>= 1            # Shift to get the next bit
        end
        f = MinkowskiFunctional(state)
        counts[f] += 1
    end

    DensityOfStates(n, counts)
end

Base.getindex(Ω::DensityOfStates, key::K) where {K} =
    get(Ω.data, key, 0)

Base.length(Ω::DensityOfStates) = length(Ω.data)

function Base.show(io::IO, Ω::DensityOfStates)
    print(io, "Density of States for a system size of $(Ω.n) x $(Ω.n).")
end

"""
    function save_density_of_states(dos::DensityOfStates, filename::String)

This saves a Density of States object.
"""
function save_density_of_states(dos::DensityOfStates, filename::String)
    open(filename, "w") do file
        serialize(file, dos)
    end
end

"""
    function load_density_of_states(filename::String)::DensityOfStates

This opens a Density of States object.
"""
function load_density_of_states(filename::String)::DensityOfStates
    open(filename, "r") do file
        return deserialize(file)
    end
end


struct MinkowskiDistribution
    n::Int
    λ::Int
    ρ::Int
    P::Accumulator
    function MinkowskiDistribution(Ω::DensityOfStates, λ, ρ)
        p = 1 - cdf(Distributions.Poisson(λ), ρ-1)

        distribution = Accumulator{MinkowskiFunctional, Float64}()

        for (key, value) in Ω.data
            distribution[key] += value * p^key.A * (1 - p)^(Ω.n^2 - key.A)
        end

        new(Ω.n, λ, ρ, distribution)
    end
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

function deviation_strength(h0_distribution::MinkowskiDistribution, x::MinkowskiFunctional)
    p = pdf(h0_distribution, x)
    ps = pdf(h0_distribution)
    mask = ps .<= p
    return -log10(sum(ps[mask]))
end

function marginalize(P::MinkowskiDistribution, field)
    marginalized_distribution = Accumulator{Int64, Float64}()

    for (x, p) in P.P
        marginalized_distribution[getfield(x, field)] += p
    end


    DiscreteNonParametric(collect(keys(marginalized_distribution)), collect(values(marginalized_distribution)))
end
