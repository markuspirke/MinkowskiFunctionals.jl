using Test
using MinkowskiFunctionals

using Distributions

@testset "generate" begin
    n = 2
    λ = 10
    ρ = 8
    d_A, d_P, d_χ = generate_distributions(n, λ, ρ)
    @test sum(probs(d_P)) ≈ 1
    @test sum(probs(d_χ)) ≈ 1
    odd_probs = probs(d_P)[isodd.(support(d_P))]
    @test odd_probs == zeros(length(odd_probs)) # perimeter cannot be odd

    n = 3
    d_A, d_P, d_χ = generate_distributions(n, λ, ρ)
    @test sum(probs(d_P)) ≈ 1
    @test sum(probs(d_χ)) ≈ 1
    odd_probs = probs(d_P)[isodd.(support(d_P))]
    @test odd_probs == zeros(length(odd_probs)) # perimeter cannot be odd
end
