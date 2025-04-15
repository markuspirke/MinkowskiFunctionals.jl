using Test
using MinkowskiFunctionals

@testset "types" begin
    counts_map = CountsMap([1 2; 3 4])
    @test 1 == counts_map[1, 1]
    @test (2, 2) == size(counts_map)
    @test (2, 3) == size(CountsMap((2, 3), 1.0))
    @test (2, 2) == size(CountsMap(2, 1.0))
    counts_map2 = CountsMap([2 4; 6 8])
    y = counts_map + counts_map
    @test counts_map2.pixels == y.pixels


    int_map = IntensityMap([1.0 2.0; 3.0 4.0])
    @test 1.0 == int_map[1, 1]
    @test (2, 2) == size(int_map)

    counts_map = CountsMap(int_map)
    @test (2, 2) == size(counts_map)

    counts_map = CountsMap([1 2; 3 4])
    bw_map = BWMap(counts_map, 3)
    @test 3 == bw_map.ρ
    @test [0 0; 1 1] == bw_map.pixels

end
