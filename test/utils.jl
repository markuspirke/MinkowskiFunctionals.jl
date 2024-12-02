using Test
using MinkowskiFunctionals
using Distributions

@testset "utils" begin
    gaussian = Normal(0, 1)
    one_σ = 2*(cdf(gaussian, 1.0) - 0.5)
    @test p2σ(1 - one_σ) ≈ 1
    five_σ = 2*(cdf(gaussian, 5.0) - 0.5)
    @test p2σ(1 - five_σ) ≈ 5

    five_σ = 2*(cdf(gaussian, 5.0) - 0.5)
    one_σ = 2*(cdf(gaussian, 1.0) - 0.5)

    @test σ2p(5.0) ≈ 1 - five_σ
    @test σ2p(1.0) ≈ 1 - one_σ


end
