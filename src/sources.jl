"""
    function flat_background(n::Int, λ::Float64)

Generates an `n × n` `CountsMap` sampled from a homogeneous Poisson field
with rate `λ`. To simulate a slightly over-background field use e.g. `1.1 * λ`.
"""
function flat_background(n::Int, λ::Float64)
    return CountsMap(n, λ)
end

"""
    function gaussian_source(n::Int, λ::Float64, A::Float64, σ::Float64)

Generates an `n × n` `CountsMap` with a homogeneous Poisson background of rate `λ`
and a Gaussian source centred at the middle of the map. `A` is the peak amplitude
(extra expected counts at the centre pixel) and `σ` is the width in pixels.
Each pixel is sampled independently from `Poisson(λ + A · exp(−r²/(2σ²)))`.
"""
function gaussian_source(n::Int, λ::Float64, A::Float64, σ::Float64)
    cx = (n + 1) / 2
    cy = (n + 1) / 2
    rates = [λ + A * exp(-((i - cx)^2 + (j - cy)^2) / (2σ^2)) for i in 1:n, j in 1:n]
    return CountsMap(IntensityMap(rates))
end

"""
    function point_source(n::Int, λ::Float64, S::Float64)

Generates an `n × n` `CountsMap` with a homogeneous Poisson background of rate `λ`
and a point source of strength `S` at the centre pixel. All pixels are sampled from
`Poisson(λ)` except the centre, which is sampled from `Poisson(λ + S)`.
"""
function point_source(n::Int, λ::Float64, S::Float64)
    rates = fill(λ, n, n)
    cx = (n + 1) ÷ 2
    cy = (n + 1) ÷ 2
    rates[cx, cy] += S
    return CountsMap(IntensityMap(rates))
end

"""
    function gradient_source(n::Int, λ::Float64, α::Float64)

Generates an `n × n` `CountsMap` with a linear intensity gradient across the map.
The rate varies as `λ · (1 + α · (2j/(n-1) - 1))` along the column direction,
so the mean rate is exactly `λ` and the total expected counts equal `n² · λ` —
indistinguishable from flat background by Li-Ma, but detectable by Minkowski analysis.
`α ∈ [0, 1)` controls the gradient strength: 0 is flat, approaching 1 gives rates
ranging from near 0 on one side to `2λ` on the other.
"""
function gradient_source(n::Int, λ::Float64, α::Float64)
    rates = [λ * (1.0 + α * (2(j - 1) / (n - 1) - 1.0)) for i in 1:n, j in 1:n]
    return CountsMap(IntensityMap(rates))
end
