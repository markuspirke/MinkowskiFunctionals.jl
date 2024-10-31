functionals = [
    (0.0, 0.0, 0.0), (1/4, 1.0, 1/4), (1/4, 1.0, 1/4), (1/2, 1.0, 0.0),
    (1/4, 1.0, 1/4), (1/2, 1.0, 0.0), (1/2, 2.0, -1/2), (3/4, 1.0, -1/4),
    (1/4, 1.0, 1/4), (1/2, 2.0, -1/2), (1/2, 1.0, 0.0), (3/4, 1.0, -1/4),
    (1/2, 1.0, 0.0), (3/4, 1.0, -1/4), (3/4, 1.0, -1/4), (1.0, 0.0, 0.0)
]
bases = [SMatrix{2, 2}(digits(UInt8(i), base=2, pad=2^2)) for i in 0:15]

functionals_old = [MinkowskiFunctional(f...) for f in functionals]
bases_states = [BasesState(b, f) for (b, f) in zip(bases, functionals_old)]


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
                if y[i, j] == bs[k][1, 1] && y[i, j+1] == bs[k][1, 2] && y[i+1, j] == bs[k][2, 1] && y[i+1, j+1] == bs[k][2, 2]
                    A += fs[k][1]
                    P += fs[k][2]
                    χ += fs[k][3]
                end
            end
        end
    end

    return MinkowskiFunctional(A, P, χ) # +1
end

MinkowskiFunctional(x::BWMap) = MinkowskiFunctional(x.pixels)

function MinkowskiFunctional_old(x::T, bases_states=bases_states) where T <: AbstractArray
    A = 0 # +1
    P = 0 # +1
    χ = 0 # +1
    y = PaddedView(0, x, size(x) .+2, (2,2))
    for j in 1:size(y)[2]-1
        for i in 1:size(y)[1]-1
            xview = @view y[i:i+1, j:j+1]
            for b in bases_states
                if xview == b.pixels
                    A += b.functional.A
                    P += b.functional.P
                    χ += b.functional.χ
                end
            end
        end
    end

    return MinkowskiFunctional(A, P, χ) # +1
end

