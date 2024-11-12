module MinkowskiFunctionals

using Base.Threads
import Base: getindex
using Distributions
using PaddedViews
using Random
using StaticArrays
using StatsBase


export CountsMap, BWMap, MinkowskiFunctional,
    PoissonMinkowskiDistributions,
    rand, Bernoulli, sample_functionals, SampledPoissonMinkowskiDistributions #,
    # MinkowskiMap, deviation_strength


include("types.jl")
include("calculate.jl")
include("generate.jl")
include("sampling.jl")
# include("minkwoski_map.jl")


end
