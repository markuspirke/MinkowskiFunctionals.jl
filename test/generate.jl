using Test
using MinkowskiFunctionals

@testset "generate" begin
    n = 2
    λ = 10
    ρ = 8
    d_A, d_P, d_χ = generate_distributions(n, λ, ρ)
    @test sum(values(d_P)) ≈ 1
    @test sum(values(d_χ)) ≈ 1
    odd_P = [d_P[i] for i in 1:2:4n^2]
    @test odd_P == [0.0 for _ in 1:2:4n^2] # perimeter cannot be odd

    n = 3
    d_A, d_P, d_χ = generate_distributions(3, λ, ρ)
    @test sum(values(d_P)) ≈ 1
    @test sum(values(d_χ)) ≈ 1
    odd_P = [d_P[i] for i in 1:2:4n^2]
    @test odd_P == [0.0 for _ in 1:2:4n^2] # perimeter cannot be odd
end
