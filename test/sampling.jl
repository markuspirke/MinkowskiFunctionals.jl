using Test
using MinkowskiFunctionals

@testset "sampling" begin
    N = 100000
    n = 2
    λ = 10
    ρ = 10
    sampled_distributions =  SampledPoissonMinkowskiDistributions(N, n, λ, ρ)
    correct_distributions = PoissonMinkowskiDistributions(n, λ, ρ)

    # not a perfect test but at least something
    @test round(mean(sampled_distributions.A), digits=1) ≈ round(mean(correct_distributions.A), digits=1)
    @test round(std(sampled_distributions.A), digits=1) ≈ round(std(correct_distributions.A), digits=1)

    @test 1.0 ≈ sum(probs(sampled_distributions.A))
    @test 1.0 ≈ sum(probs(sampled_distributions.P))
    @test 1.0 ≈ sum(probs(sampled_distributions.χ))
end
