struct CountsMap
    pixels
end

function CountsMap(d::FullNormal, n, N)
    samples = rand(d, N)
    counts_map = fit(Histogram, (samples[1, :], samples[2,:]), (range(-1, 1, n+1), range(-1,1,n+1)))

    return CountsMap(counts_map.weights)
end

struct BWMap
    ρ
    pixels::Matrix{Bool}
end


function BWMap(x::CountsMap, ρ)
    return BWMap(ρ, x.pixels .>= ρ)
end

"""
    struct MinkowskiFunctional

This is a data type for MinkowskiFunctionals.
"""
struct MinkowskiFunctional
    A
    P
    χ
end

function Base.show(io::IO, m::MinkowskiFunctional)
    print(io, "A: $(m.A), P: $(m.P), χ: $(m.χ)")
end
