"""
    function generate_distributions(n, λ, ρ)

This generates for a given n, λ and ρ the distributions of Minkowski functionals.
The distribution for the area A is returned as a Binomial distribution from the Distributions.jl pacakge.
The distributions for perimeter P and Euler characteristic are returned as dictionaries.
Where the keys are P and χ and the values are the corresponding probabilities.
"""
function generate_distributions(n, λ, ρ)
    n_microstates = 2^(n^2)
    d_poisson = Distributions.Poisson(λ)
    p = 1 - cdf(d_poisson, ρ-1)

    # lower and upper bounds for the Minkowski functionals in 2D
    l_area = 0
    u_area = n^2
    l_perimeter = 0
    u_perimeter = 4n^2
    l_euler = 2 - n^2
    u_euler = isodd(n) ? ((n+1)/2)^2 : (n/2)^2
    N_euler = u_euler - l_euler


    d_area = Binomial(n^2, p)
    d_perimeter = Dict(zip(l_perimeter:u_perimeter, [0.0 for _ in 1:(u_perimeter - l_perimeter)]))
    d_euler = Dict(zip(l_euler:u_euler, [0.0 for _ in 1:(u_euler - l_euler + 1)]))


    state = BitMatrix(undef, n, n)

    # Iterate over all 2^(n^2) binary configurations
    for i in 0:2^(n^2) - 1
        # Fill the matrix `state` using bitwise operations
        x = i
        for k in 1:n^2
            state[k] = x & 1  # Set each element based on the binary representation
            x >>= 1            # Shift to get the next bit
        end
        f = MinkowskiFunctional(state)
        p_f = pdf(d_area, f.A)
        normalization = binomial(n^2, Int(f.A))
        d_perimeter[f.P] += p_f/normalization
        d_euler[f.χ] += p_f/normalization
    end

    return d_area, d_perimeter, d_euler
end

function foo(n)
    functionals = MinkowskiFunctional[]
    for i in 1:2^(n^2)
        state = SMatrix{n, n}(digits(UInt32(i-1), base=2, pad=n^2))
        push!(functionals, MinkowskiFunctional(state))
    end

    functionals
end

function foo_optimized(n)
    # Pre-allocate a mutable n x n matrix that will be reused
    # states = [BitMatrix(undef, n, n) for _ in 1:2^(n^2)]
    functionals = MinkowskiFunctional[]
    state = BitMatrix(undef, n, n)

    # Iterate over all 2^(n^2) binary configurations
    for i in 0:2^(n^2) - 1
        # Fill the matrix `state` using bitwise operations
        x = i
        for k in 1:n^2
            state[k] = x & 1  # Set each element based on the binary representation
            x >>= 1            # Shift to get the next bit
        end
        # `state` now contains the matrix in its current configuration.
        # Do something with `state` if needed, like computation or printing.
        push!(functionals, MinkowskiFunctional(state))
        # push!(states, state)
    end

    functionals
end
