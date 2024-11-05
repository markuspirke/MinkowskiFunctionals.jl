# TABLE FROM MOPRHOMETRIC ANALYSIS PAPER 2013
functionals = [
    (0.0, 0.0, 0.0), (1/4, 1.0, 1/4), (1/4, 1.0, 1/4), (1/2, 1.0, 0.0),
    (1/4, 1.0, 1/4), (1/2, 1.0, 0.0), (1/2, 2.0, -1/2), (3/4, 1.0, -1/4),
    (1/4, 1.0, 1/4), (1/2, 2.0, -1/2), (1/2, 1.0, 0.0), (3/4, 1.0, -1/4),
    (1/2, 1.0, 0.0), (3/4, 1.0, -1/4), (3/4, 1.0, -1/4), (1.0, 0.0, 0.0)
]
bases = [SMatrix{2, 2}(digits(UInt8(i), base=2, pad=2^2)) for i in 0:15]

"""
    function MinkowskiFunctional(x::T, bs=bases, fs=functionals) where T <: AbstractArray

For a given black and white image x, this calculates the Minkowski functionals
for the whole image based on a 2 x 2 look up table.
"""
function MinkowskiFunctional(x::T, bs=bases, fs=functionals) where T <: AbstractArray
    A = 0.0
    P = 0.0
    χ = 0.0

    y = PaddedView(0, x, size(x) .+2, (2,2))
    m, n = size(y)
    for j in 1:n-1
        for i in 1:m-1
            for k in 1:length(bs)
                @inbounds if y[i, j] == bs[k][1, 1] && y[i, j+1] == bs[k][1, 2] && y[i+1, j] == bs[k][2, 1] && y[i+1, j+1] == bs[k][2, 2]
                    @inbounds A += fs[k][1]
                    @inbounds P += fs[k][2]
                    @inbounds χ += fs[k][3]
                    break
                end
            end
        end
    end

   return MinkowskiFunctional(A, P, χ)
end

MinkowskiFunctional(x::BWMap) = MinkowskiFunctional(x.pixels)

