module MinkowskiFunctionals

using Base.Threads
import Base: getindex, +, append!, size, show
using Distributions
using DataStructures
using PaddedViews
using Random
using StaticArrays
using StatsBase
using HDF5
using ProgressMeter
using Dates


export CountsMap, IntensityMap, BWMap, MinkowskiFunctional, IntX, p2σ, σ2p, Background,
    PoissonMinkowskiDistributions, DensityOfStates,
    MinkowskiDistribution, pdf, marginalize, reduce_functional, AreaDistribution, append!,
    add_functionals!, compatibility,
    save_density_of_states, load_density_of_states, convert_counter,
    rand, Bernoulli, sample_functionals, SampledPoissonMinkowskiDistributions,
    MinkowskiMap, minkowski_map_A, correct_trials, get_tresholds,
    λ_lima, significance_lima, lima_map, lima_map_roundkernel


include("types.jl")
include("utils.jl")
include("calculate.jl")
include("generate.jl")
include("density_of_states.jl")
include("distributions.jl")
include("sampling.jl")
include("minkowski_map.jl")
include("lima.jl")


end
