"""
This generates a Minkowski map from a given CountsMap.
This also needs distributions of Minkowski functionals
and one needs to decide which functional is used.
"""
struct MinkowskiMap
    λ
    ρ
    L
    pixels

    function MinkowskiMap(x, h0_distributions, F)
        m, n = size(x.pixels)
        ρs = [get_distribution(s).ρ for s in h0_distributions]
        λ = get_distribution(h0_distributions[1]).λ
        L = get_distribution(h0_distributions[1]).n
        l = floor(Int, L/2)
        Ds = zeros(n - 2l, m - 2l)
        for j in l+1:m-l
            for i in l+1:n-l
                deviation_strengths = zeros(length(ρs))
                for (k, ρ) in enumerate(ρs)
                    bw_map = BWMap(x, ρ)
                    functional = MinkowskiFunctional(bw_map.pixels[i-l:i+l, j-l:j+l])
                    deviation_strengths[k] = deviation_strength(getfield(get_distribution(h0_distributions[k]), F), getfield(functional, F))
                end
                Ds[i-l, j-l] = maximum(deviation_strengths)
            end
        end

        return new(λ, ρs, L, Ds)
    end
end


"""
    function deviation_strength(h0_distribution::T, x) where {T<:DiscreteDistribution}

This takes a discrete distribution and calulates the deviation strength for
the event X, where X needs to be in the support (sample space) of the distribution.
"""
function deviation_strength(h0_distribution::T, x) where {T<:DiscreteDistribution}
    p = pdf(h0_distribution, x)
    # if p == 0
    #     @show p, x
    # end
    mask = probs(h0_distribution) .<= p
    return -log10(sum(probs(h0_distribution)[mask]))
end