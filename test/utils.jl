using Test
using MinkowskiFunctionals
using Distributions

@testset "utils" begin
    gaussian = Normal(0, 1)
    five_sigma = 1 - 2*(cdf(gaussian, 5.0) - 0.5)
    @test deviation2σ(-log10(five_sigma)) ≈ 5
end
