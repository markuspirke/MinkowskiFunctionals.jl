module MinkowskiFunctionals

using Base.Threads
import Base: getindex, +
using Distributions
using DataStructures
using Serialization
using PaddedViews
using Random
using StaticArrays
using StatsBase


export CountsMap, BWMap, MinkowskiFunctional, IntX,
    PoissonMinkowskiDistributions, DensityOfStates, MinkowskiDistribution, pdf, marginalize,
    add_functionals!,
    save_density_of_states, load_density_of_states, convert_counter,
    rand, Bernoulli, sample_functionals, SampledPoissonMinkowskiDistributions,
    MinkowskiMap, deviation_strength


include("types.jl")
include("calculate.jl")
include("generate.jl")
include("density_of_states.jl")
include("sampling.jl")
include("minkwoski_map.jl")


end