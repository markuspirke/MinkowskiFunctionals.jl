using MinkowskiFunctionals
using Test

@testset "MinkowskiFunctionals.jl" begin
    include("calculate.jl")
    include("generate.jl")
    include("density_of_states.jl")
    include("sampling.jl")
    include("minkowski_map.jl")
    # Write your tests here.
end
