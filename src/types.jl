struct CountsMap
    pixels::Matrix{Int64}
end

Base.getindex(x::CountsMap, i, j) = x.pixels[i, j]


function CountsMap(d::FullNormal, n, N)
    samples = rand(d, N)
    counts_map = fit(Histogram, (samples[1, :], samples[2,:]), (range(-1, 1, n+1), range(-1,1,n+1)))

    return CountsMap(counts_map.weights)
end

struct BWMap
    ρ::Int64
    pixels::Matrix{Bool}
end


function BWMap(x::CountsMap, ρ)
    return BWMap(ρ, x.pixels .>= ρ)
end

function BWMap(x::Matrix{Int64}, ρ)
    return BWMap(ρ, x .>= ρ)
end
"""
    struct MinkowskiFunctional

This is a data type for MinkowskiFunctionals.
"""
struct MinkowskiFunctional
    A::Int64
    P::Int64
    χ::Int64
end

function Base.show(io::IO, m::MinkowskiFunctional)
    print(io, "A: $(m.A), P: $(m.P), χ: $(m.χ)")
end


struct IntX <: Number
    base::Int64
    exp::Int64
end
Base.convert(T::Type{IntX}, x::Number) = T(x, 0)::T
Base.convert(T::Type{Int64}, x::IntX) = T(x.base) * T(T(10)^(x.exp))::T
Base.convert(T::Type{Int128}, x::IntX) = T(x.base) * T(T(10)^(x.exp))::T
Base.convert(T::Type{BigInt}, x::IntX) = T(x.base) * T(T(10)^(x.exp))::T

+(a::IntX, b::IntX) = IntX(a.base + b.base, a.exp + b.exp)
