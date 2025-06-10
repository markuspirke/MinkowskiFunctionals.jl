"""
    struct CountsMap

This is a data type that stores the counts for some region of space.
"""
struct CountsMap
    pixels::Matrix{Int64}
end

Base.size(x::CountsMap) = size(x.pixels)
Base.getindex(x::CountsMap, i, j) = x.pixels[i, j]
Base.setindex!(x::CountsMap, i, j, k) = setindex!(x.pixels, i, j, k)
+(x::CountsMap, y::CountsMap) = CountsMap(x.pixels + y.pixels)

function CountsMap(x::Tuple{Int64, Int64}, λ::Float64)
    m, n = x
    d = Poisson(λ)
    CountsMap(rand(d, m, n))
end

function CountsMap(x::Int64, λ::Float64)
    d = Poisson(λ)
    CountsMap(rand(d, x, x))
end

function CountsMap(x::Int64, d::Poisson{Float64})
    CountsMap(rand(d, x, x))
end
"""
    struct IntensityMap

This is a data type for the intensity map for a given Poisson random field.
"""
struct IntensityMap
    λs::Matrix{Float64}
end

Base.size(x::IntensityMap) = size(x.λs)
Base.getindex(x::IntensityMap, i, j) = x.λs[i, j]

function Base.show(io::IO, x::IntensityMap)
    println(io, "$(size(x)[1]) x $(size(x)[2]) Intensity Map:")
    for row in eachrow(x.λs)
        println(io, row)
    end
end

function CountsMap(x::IntensityMap)
    m, n = size(x)
    y = zeros(m, n)
    for j in 1:n
        for i in 1:m
            y[i, j] = rand(Poisson(x[i, j]))
        end
    end
    return CountsMap(y)
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

struct Background
    pixels::Matrix{Float64}
end

Base.size(x::Background) = size(x.pixels)
Base.getindex(x::Background, i, j) = x.pixels[i, j]

function CountsMap(x::Background)
    m, n = size(x)
    y = zeros(m, n)
    for j in 1:n
        for i in 1:m
            y[i, j] = rand(Poisson(x[i, j]))
        end
    end
    return CountsMap(y)
end
"""
    struct MinkowskiFunctional

This is a data type for MinkowskiFunctionals.
"""
struct MinkowskiFunctional
    A::Int16
    P::Int16
    χ::Int16
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

