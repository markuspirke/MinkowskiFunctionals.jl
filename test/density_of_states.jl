using Test
using MinkowskiFunctionals

const SAMPLES_DIR = joinpath(@__DIR__, "samples")

@testset "density_of_states" begin
    Ω = DensityOfStates(2)
    @test 1 == Ω[MinkowskiFunctional(0, 0, 0)]
    @test 4 == Ω[MinkowskiFunctional(1, 4, 1)]
    @test 4 == Ω[MinkowskiFunctional(2, 6, 1)]
    @test 2 == Ω[MinkowskiFunctional(2, 8, 1)]
    @test 4 == Ω[MinkowskiFunctional(3, 8, 1)]
    @test 1 == Ω[MinkowskiFunctional(4, 8, 1)]

    save_density_of_states(Ω, joinpath(SAMPLES_DIR, "dos"))
    Ω = load_density_of_states(joinpath(SAMPLES_DIR, "dos"))
    @test 1 == Ω[MinkowskiFunctional(0, 0, 0)]
    @test 4 == Ω[MinkowskiFunctional(1, 4, 1)]
    @test 4 == Ω[MinkowskiFunctional(2, 6, 1)]
    @test 2 == Ω[MinkowskiFunctional(2, 8, 1)]
    @test 4 == Ω[MinkowskiFunctional(3, 8, 1)]
    @test 1 == Ω[MinkowskiFunctional(4, 8, 1)]

    n, λ, ρ = 3, 10, 10
    Ω = DensityOfStates(n)
    P = MinkowskiDistribution(Ω, λ, ρ)
    @test 1.0 ≈ sum(pdf(P))

    P_A = marginalize(P, :A)
    @test 1.0 ≈ sum(pdf(P_A))
    P_P = marginalize(P, :P)
    @test 1.0 ≈ sum(pdf(P_P))
    P_χ = marginalize(P, :χ)
    @test 1.0 ≈ sum(pdf(P_χ))
end
