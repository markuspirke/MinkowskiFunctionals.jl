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

struct Edge
    a::Tuple{Int, Int}
    b::Tuple{Int, Int}
end

function find_perimeter_edges(img)
    n, m = size(img)
    edges = Edge[]
    for i in 1:n
        for j in 1:m
            if img[i, j] == 1
                if img[i-1, j] == 0
                    push!(edges, Edge((i-1, j-1), (i-1, j)))
                end
                if img[i+1, j] == 0
                    push!(edges, Edge((i, j-1), (i, j)))
                end
                if img[i, j-1] == 0
                    push!(edges, Edge((i-1, j-1), (i, j-1)))
                end
                if img[i, j+1] == 0
                    push!(edges, Edge((i-1, j), (i, j)))
                end
            end
        end
    end
    edges
end


function label_outer_background(img)
    n, m = size(img)
    visited = falses(n, m)
    queue = Tuple{Int,Int}[]

    # Start from boundary white pixels
    for i in 1:n
        for j in 1:m
            if (i in (1, n) || j in (1, m)) && img[i, j] == 0
                push!(queue, (i, j))
            end
        end
    end

    while !isempty(queue)
        (ci, cj) = popfirst!(queue)
        if ci < 1 || ci > n || cj < 1 || cj > m
            continue
        end
        if visited[ci, cj] || img[ci, cj] != 0
            continue
        end
        visited[ci, cj] = true
        append!(queue, [
            (ci-1, cj), (ci+1, cj),
            (ci, cj-1), (ci, cj+1)
        ])
    end

    visited .== 0  # true = outer background
end

function get_round_kernel(L)
    l = floor(L/2)
    cx, cy = ((L+1)/2, (L+1)/2)
    [ (i - cx)^2 + (j - cy)^2 ≤ l^2 for i in 1:L, j in 1:L ]
end
