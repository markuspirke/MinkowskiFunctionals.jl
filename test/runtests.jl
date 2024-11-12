using MinkowskiFunctionals
using Test

@testset "MinkowskiFunctionals.jl" begin
    include("calculate.jl")
    include("generate.jl")
    include("sampling.jl")
    # Write your tests here.
end
