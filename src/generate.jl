"""
Datatype for holding distributions of Minkowski functionals for nxn images.
"""
struct PoissonMinkowskiDistributions
    n
    λ
    ρ
    A
    P
    χ
end

function get_distribution(x::PoissonMinkowskiDistributions)
    return x
end
"""
    function PoissonMinkowskiDistributions(n, λ, ρ)

This generates for a given n, λ and ρ the distributions of Minkowski functionals.
The distribution for the area A is returned as a Binomial distribution from the Distributions.jl pacakge.
The distributions for perimeter P and Euler characteristic are returned as a type DiscreteParametric from
Distribution.jl, which is essentially a probability mass function.
"""
function PoissonMinkowskiDistributions(n, λ, ρ)
    n_microstates = 2^(n^2)
    d_poisson = Distributions.Poisson(λ)
    p = 1 - cdf(d_poisson, ρ-1)

    # lower and upper bounds for the Minkowski functionals in 2D
    l_perimeter = 0
    u_perimeter = 4n^2
    l_euler = 2 - n^2
    u_euler = isodd(n) ? ((n+1)/2)^2 : (n/2)^2


    d_area = Binomial(n^2, p)
    dict_perimeter = Dict(zip(l_perimeter:u_perimeter, [0.0 for _ in 1:(u_perimeter - l_perimeter)]))
    dict_euler = Dict(zip(l_euler:u_euler, [0.0 for _ in 1:(u_euler - l_euler + 1)]))


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
        dict_perimeter[f.P] += p_f/normalization
        dict_euler[f.χ] += p_f/normalization
    end

    d_perimeter = DiscreteNonParametric(collect(keys(dict_perimeter)), collect(values(dict_perimeter)))
    d_euler = DiscreteNonParametric(collect(keys(dict_euler)), collect(values(dict_euler)))


    return PoissonMinkowskiDistributions(n, λ, ρ, d_area, d_perimeter, d_euler)
end





