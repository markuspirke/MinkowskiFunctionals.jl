using MinkowskiFunctionals
using Test

@testset "MinkowskiFunctionals.jl" begin
    include("types.jl")
    include("calculate.jl")
    include("generate.jl")
    include("density_of_states.jl")
    include("distributions.jl")
    include("sampling.jl")
    include("minkowski_map.jl")
    include("utils.jl")
    include("lima.jl")
    include("test_statistic.jl")
    # Write your tests here.
end
