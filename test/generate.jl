using Test
using MinkowskiFunctionals

using Distributions

@testset "generate" begin
    n = 2
    λ = 10
    ρ = 8
    d = PoissonMinkowskiDistributions(n, λ, ρ)
    @test sum(pdf(d.A)) ≈ 1
    @test sum(pdf(d.P)) ≈ 1
    @test sum(pdf(d.χ)) ≈ 1
    odd_probs = pdf(d.P)[isodd.(support(d.P))]
    @test odd_probs == zeros(length(odd_probs)) # perimeter cannot be odd

    n = 3
    d = PoissonMinkowskiDistributions(n, λ, ρ)
    @test sum(pdf(d.A)) ≈ 1
    @test sum(pdf(d.P)) ≈ 1
    @test sum(pdf(d.χ)) ≈ 1
    odd_probs = pdf(d.P)[isodd.(support(d.P))]
    @test odd_probs == zeros(length(odd_probs)) # perimeter cannot be odd
end
