"""
    function p2σ(D)

Converts a probability to a significance as given in σ of a unit Gaussian.
"""
function p2σ(p)
	if p > 1e-15
		q = 1 - p
		return quantile(Normal(0, 1), 1/2+q/2) # till p=0.5 we have 0, then go for half of p again
	else
        return sqrt(-2*log(p))
    end
end

"""
    function σ2p(D)

Converts a a significance as given in σ of a unit Gaussian to a probability.
"""
function σ2p(σ)
	1 - 2 * (cdf(Normal(0, 1), σ) - 0.5)
end

"""
    function binary_search(xs::Vector{T}, x::T) where T

A self implemented binary search.
"""
function binary_search(xs::Vector{T}, x::T) where T
    left = 1
    right = length(xs)
    mid = 1

    while left <= right
        mid = (left + right) ÷ 2
        @inbounds if xs[mid] ≈ x
            return mid
        elseif xs[mid] < x
            left = mid + 1
        else
            right = mid - 1
        end
    end

    return mid
end
