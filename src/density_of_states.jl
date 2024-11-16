struct DensityOfStates
    n
    data::Accumulator{MinkowskiFunctional, Int64}
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
