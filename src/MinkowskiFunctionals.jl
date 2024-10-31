module MinkowskiFunctionals

# Write your package code here.
using Distributions
using StatsBase
using StaticArrays
using PaddedViews
import Base: rand
export CountsMap, BWMap, MinkowskiFunctional
include("types.jl")
include("calculate.jl")

end
