"""
    function Base.rand(T::Type{CountsMap}, d::Poisson, n)
using StatsBase: sample_ordered!

New method for the rand function to generate random Poisson counts maps.
"""
function Random.rand(T::Type{CountsMap}, d::Poisson, n)
    return T(rand(d, n, n))
end


"""
    function Base.rand(T::Type{BWMap}, d::Binomial, n)

New method for the rand function to generate random black and white maps for a given
threshold and a given Poisson background.
"""
function Random.rand(T::Type{BWMap}, d::Bernoulli, n) # THIS IS BAD
    return T(rand(d, n, n)) # ONE WOULD EXPECT TO GET A SOMETHING DIFFERENT ACCORDING TO RANDOM I THINK
end

function _random_bwmap(d::Bernoulli, n)
    return rand(d, n, n)
end
"""
    function Distributions.Binomial(n, λ, ρ)

Helper function which generates the correct Binomial distribution needed for sampling
black and white maps.
Here,
 - n is the dimension of the image
 - λ the expected background rate
 - ρ the threshold that is chosen
"""
function Distributions.Bernoulli(λ, ρ)
    d = Poisson(λ)
    p = 1 - cdf(d, ρ-1)

    return Bernoulli(p)
end

"""
    function sample_functionals(N, λ, ρ, n)

This samples N nxn black and white images for a given λ and ρ and calculates the Minkowski functionals.
"""
function sample_functionals(N, n, λ, ρ)
    d = Bernoulli(λ, ρ)
    functionals = Vector{MinkowskiFunctional}(undef, N)
    Threads.@threads for i in 1:N
        functionals[i] = MinkowskiFunctional(_random_bwmap(d, n))
    end

    functionals
end

struct SampledPoissonMinkowskiDistributions
    N
    d::PoissonMinkowskiDistributions

    function SampledPoissonMinkowskiDistributions(N, n, λ, ρ)
        functionals = sample_functionals(N, n, λ, ρ)
        new(N, _minkowski_distributions(functionals, n, λ, ρ))
    end
end

function _minkowski_distributions(functionals::Vector{MinkowskiFunctional}, n, λ, ρ)
    l_area = 0
    u_area = n^2
    l_perimeter = 0
    u_perimeter = 4n^2
    l_euler = 2 - n^2
    u_euler = isodd(n) ? ((n+1)/2)^2 : (n/2)^2


    dict_area = Dict(zip(l_area:u_area, [0.0 for _ in 1:(u_area - l_area + 1)]))
    dict_perimeter = Dict(zip(l_perimeter:u_perimeter, [0.0 for _ in 1:(u_perimeter - l_perimeter)]))
    dict_euler = Dict(zip(l_euler:u_euler, [0.0 for _ in 1:(u_euler - l_euler + 1)]))

    for functional in functionals
        dict_area[functional.A] += 1
        dict_perimeter[functional.P] += 1
        dict_euler[functional.χ] += 1
    end

    norm = length(functionals)
    d_area = DiscreteNonParametric(collect(keys(dict_area)), collect(values(dict_area))./norm)
    d_perimeter = DiscreteNonParametric(collect(keys(dict_perimeter)), collect(values(dict_perimeter))./norm)
    d_euler = DiscreteNonParametric(collect(keys(dict_euler)), collect(values(dict_euler))./norm)

    return PoissonMinkowskiDistributions(n, λ, ρ, d_area, d_perimeter, d_euler)
end

function get_distribution(x::SampledPoissonMinkowskiDistributions)
    return x.d
end
