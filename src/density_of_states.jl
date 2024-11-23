struct DensityOfStates
    n::Int
    data::Accumulator{MinkowskiFunctional, Number}
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

"""
    function add_functionals!(counter::Accumulator{MinkowskiFunctional, T}, fname::AbstractString) where {T<:Number}

This reads a density of states .dat file and updates a given counter.
"""
function add_functionals!(counter::Accumulator{MinkowskiFunctional, T}, fname::AbstractString) where {T<:Number}
    open(fname, "r") do f
        A = split(split(fname, "_A_")[end], "_")[1]
        A = parse(Int64, A)
        for line in eachline(f)
            x = replace(line, r"  +" => ' ')
            if occursin("e", x)
                x = split(x, " ")
                P, χ, N = x[end-2], x[end-1], x[end]
                P, χ, N = parse(T, P), parse(T, χ), parse(Float64, N)
                N = round(T, N)
                functional = MinkowskiFunctional(A, P, χ)
                counter[functional] += N
            else
                x = split(x, " ")
                P, χ, N = x[end-2], x[end-1], x[end]
                P, χ, N = parse(T, P), parse(T, χ), parse(T, N)
                functional = MinkowskiFunctional(A, P, χ)
                counter[functional] += N
            end
        end
    end
end

"""
    function add_functionals!(counter::Accumulator{MinkowskiFunctional, IntX}, fname::AbstractString)

This reads a density of states .dat file and updates a given counter.
"""
function add_functionals!(counter::Accumulator{MinkowskiFunctional, IntX}, fname::AbstractString)
    open(fname, "r") do f
        A = split(split(fname, "_A_")[end], "_")[1]
        A = parse(Int64, A)
        for line in eachline(f)
            x = replace(line, r"  +" => ' ')
            x = split(x, " ")
            P, χ, N = x[end-2], x[end-1], x[end]
            P, χ = parse(Int, P), parse(Int, χ)
            if occursin("e+", N)
                base, exp = split(N, "e+")
                exp = parse(Int, exp)
                if length(base)-2 < exp
                    exp = exp - (length(base) - 2)
                    base = parse(Int, replace(base, r"\." => ""))
                    N = IntX(base, exp)
                else
                    base = round(parse(Float64, base), digits=exp)
                    base = floor(Int, base*10^exp)
                    exp = 0
                    N = IntX(base, exp)
                end
            elseif occursin("e-", N)
                base = Int(round(parse(Float64, N), digits=1))
                exp = 0
                N = IntX(base, exp)
            else
                base = parse(Int, N)
                exp = 0
                N = IntX(base, exp)
            end
            functional = MinkowskiFunctional(A, P, χ)
            counter[functional] += N
        end
    end
end

"""
    function convert_counter(T::Type{S}, counter::Accumulator{MinkowskiFunctional, IntX}) where {S<:Number}

This converts a Counter of IntX to the best Int counter.
"""
function convert_counter(T::Type{S}, counter::Accumulator{MinkowskiFunctional, IntX}) where {S<:Number}
    new_counter = Accumulator{MinkowskiFunctional, T}()
    for (key, value) in counter
        new_counter[key] += convert(T, value)
    end
    return new_counter
end
"""
    function DensityOfStates(dirname::AbstractString)

This takes the path to one of the directories inside the parameters directory.
Then the density of states is calculated from the files inside the directory.
"""
function DensityOfStates(dirname::AbstractString)
    pattern = r"(.)x(.)"
    n = parse(Int, match(pattern, dirname)[1])
    fnames = readdir(dirname, join=true)
    counterX = Accumulator{MinkowskiFunctional, IntX}()
    for fname in fnames
        add_functionals!(counterX, fname)
    end
    max_exp = 16 * maximum(getfield.(collect(values(counterX)), :exp))
    if max_exp < length(digits(typemax(Int64)))
        return DensityOfStates(n, convert_counter(Int64, counterX))
    elseif max_exp < length(digits(typemax(Int128)))
        return DensityOfStates(n, convert_counter(Int128, counterX))
    else
        return DensityOfStates(n, convert_counter(BigInt, counterX))
    end
end
