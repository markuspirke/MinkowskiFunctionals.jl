module MinkowskiFunctionals

using Base.Threads
import Base: getindex, +, append!, size, show, minimum, maximum, abs, download
import Statistics: mean
using Distributions
using DataStructures
using PaddedViews
using Random
using StaticArrays
using StatsBase
using HDF5
using ProgressMeter
using Dates
using SpecialFunctions

export CountsMap, IntensityMap, BWMap, MinkowskiFunctional, IntX, p2σ, σ2p, Background,
    PoissonMinkowskiDistributions, DensityOfStates,
    MinkowskiDistribution, pdf, marginalize, reduce_functional, AreaDistribution, append!, window_size,
    write_pvalues, write_necessary_pvalues, read_pvalues,
    add_functionals!, compatibility,
    save_density_of_states, load_density_of_states, convert_counter,
    rand, Bernoulli, sample_functionals, SampledPoissonMinkowskiDistributions,
    MinkowskiMap, minkowski_map_A, correct_trials, get_thresholds, get_λs,
    λ_lima, significance_lima, lima_map, lima_map_roundkernel,
    ECCDF, write_eccdf, read_eccdf, calc_ts!, calc_ts, update_distributions!, find_max_threshold


include("types.jl")
include("utils.jl")
include("calculate.jl")
include("generate.jl")
include("density_of_states.jl")
include("distributions.jl")
include("sampling.jl")
include("minkowski_map.jl")
include("lima.jl")
include("test_statistic.jl")


end
