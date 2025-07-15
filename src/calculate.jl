# TABLE FROM MOPRHOMETRIC ANALYSIS PAPER 2013
const functionals = SArray{Tuple{2, 2, 2, 2}, Tuple{Float64, Float64, Float64}, 4}([
    (0.0, 0.0, 0.0), (1/4, 1.0, 1/4), (1/4, 1.0, 1/4), (1/2, 1.0, 0.0),
    (1/4, 1.0, 1/4), (1/2, 1.0, 0.0), (1/2, 2.0, -1/2), (3/4, 1.0, -1/4),
    (1/4, 1.0, 1/4), (1/2, 2.0, -1/2), (1/2, 1.0, 0.0), (3/4, 1.0, -1/4),
    (1/2, 1.0, 0.0), (3/4, 1.0, -1/4), (3/4, 1.0, -1/4), (1.0, 0.0, 0.0)
])
"""
    function MinkowskiFunctional(x::T, fs=functionals) where T <: AbstractArray

For a given black and white image x, this calculates the Minkowski functionals
for the whole image based on a 2 x 2 look up table.
"""
function MinkowskiFunctional(x::T, fs=functionals) where T <: AbstractArray
    A = 0.0
    P = 0.0
    χ = 0.0

    y = PaddedView(0, x, size(x) .+2, (2,2))
    m, n = size(y)
    for j in 1:n-1
        for i in 1:m-1
            _A, _P, _χ = fs[y[i,j]+1, y[i,j+1]+1, y[i+1,j]+1, y[i+1,j+1]+1]  # 1-indexed arrays
            A += _A
            P += _P
            χ += _χ
        end
    end

   return MinkowskiFunctional(Int(A), Int(P), Int(χ))
end

MinkowskiFunctional(x::BWMap) = MinkowskiFunctional(x.pixels)

