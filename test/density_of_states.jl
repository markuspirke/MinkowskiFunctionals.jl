using Test
using MinkowskiFunctionals
using Distributions
using DataStructures
using HDF5

const SAMPLES_DIR = joinpath(@__DIR__, "samples")

@testset "density_of_states" begin
    n = 2
    Ω = DensityOfStates(n)
    @test 1 == Ω[MinkowskiFunctional(0, 0, 0)]
    @test 4 == Ω[MinkowskiFunctional(1, 4, 1)]
    @test 4 == Ω[MinkowskiFunctional(2, 6, 1)]
    @test 2 == Ω[MinkowskiFunctional(2, 8, 1)]
    @test 4 == Ω[MinkowskiFunctional(3, 8, 1)]
    @test 1 == Ω[MinkowskiFunctional(4, 8, 1)]
    @test 2^(n^2) == sum(values(Ω.data))

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
    @test 2^(n^2) == sum(values(Ω.data))

    # TO READ MICHAEL KLATT PRECALCULATED DOS

    counter = Accumulator{MinkowskiFunctional, Int64}()
    add_functionals!(counter, joinpath(SAMPLES_DIR, "DoS_averaged_A_2_PC_pix_wbc_5x5_seed_0_seedwght_0.dat"))

    @test counter[MinkowskiFunctional(2, 6, 1)] == 40
    @test counter[MinkowskiFunctional(2, 8, 1)] == 32
    @test counter[MinkowskiFunctional(2, 8, 2)] == 228

    counter_X = Accumulator{MinkowskiFunctional, IntX}()
    add_functionals!(counter_X, joinpath(SAMPLES_DIR, "DoS_averaged_A_2_PC_pix_wbc_5x5_seed_0_seedwght_0.dat"))

    @test counter_X[MinkowskiFunctional(2, 6, 1)] == IntX(40, 0)
    @test counter_X[MinkowskiFunctional(2, 8, 1)] == IntX(32, 0)
    @test counter_X[MinkowskiFunctional(2, 8, 2)] == IntX(228, 0)

    new_counter = convert_counter(Int64, counter_X)
    @test new_counter == counter

    counter_X = Accumulator{MinkowskiFunctional, IntX}()
    add_functionals!(counter_X, joinpath(SAMPLES_DIR, "DoS_averaged_A_46_PC_pix_wbc_11x11_seed_0_seedwght_0.dat"))
    @test counter_X[MinkowskiFunctional(46, 156, -14)] == IntX(26775834197587961, 5)
    counter = convert_counter(Int128, counter_X)
    @test counter[MinkowskiFunctional(46, 156, -14)] == 2677583419758796100000

    counter_X = Accumulator{MinkowskiFunctional, IntX}()
    add_functionals!(counter_X, joinpath(SAMPLES_DIR, "DoS_averaged_A_32_PC_pix_wbc_15x15_seed_0_seedwght_0.dat"))
    add_functionals!(counter_X, joinpath(SAMPLES_DIR, "DoS_averaged_A_170_PC_pix_wbc_15x15_seed_0_seedwght_0.dat"))
    @test counter_X[MinkowskiFunctional(32, 100, 11)] == IntX(24974287959444129, 19)
    @test counter_X[MinkowskiFunctional(170, 200, -20)] == IntX(22246633333168453, 35)
    counter = convert_counter(BigInt, counter_X)
    @test counter[MinkowskiFunctional(32, 100, 11)] == 249742879594441290000000000000000000
    @test counter[MinkowskiFunctional(170, 200, -20)] == 2224663333316845300000000000000000000000000000000000


    Ω = DensityOfStates(joinpath(SAMPLES_DIR, "structure_5x5"))
    @test Ω.n == 5
    @test Ω.data[MinkowskiFunctional(2, 6, 1)] == 40
    @test Ω.data[MinkowskiFunctional(2, 8, 1)] == 32
    @test Ω.data[MinkowskiFunctional(2, 8, 2)] == 228

    Ω = DensityOfStates(joinpath(SAMPLES_DIR, "structure_10x10"))
    @test Ω.n == 10
end
