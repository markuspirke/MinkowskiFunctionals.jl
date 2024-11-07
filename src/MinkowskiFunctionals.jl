module MinkowskiFunctionals

using Base.Threads
using Distributions
using PaddedViews
using Random
using StaticArrays
using StatsBase


export CountsMap, BWMap, MinkowskiFunctional,
    rand, Bernoulli, sample_functionals,
    generate_distributions


include("types.jl")
include("calculate.jl")
include("sampling.jl")
include("generate.jl")


end
