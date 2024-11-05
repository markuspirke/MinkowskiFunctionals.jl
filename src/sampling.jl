"""
    function Base.rand(T::Type{CountsMap}, d::Poisson, n)

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
function sample_functionals(N, λ, ρ, n)
    d = Bernoulli(λ, ρ)
    functionals = Vector{MinkowskiFunctional}(undef, N)
    Threads.@threads for i in 1:N
        functionals[i] = MinkowskiFunctional(_random_bwmap(d, n))
    end

    functionals
end
