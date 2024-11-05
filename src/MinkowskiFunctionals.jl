module MinkowskiFunctionals

using Distributions
using StatsBase
using StaticArrays
using PaddedViews
# import Random: rand
using Random
using Base.Threads


export CountsMap, BWMap, MinkowskiFunctional,
    rand, Bernoulli, sample_functionals,
    generate_distributions


include("types.jl")
include("calculate.jl")
include("sampling.jl")
include("generate.jl")


end
