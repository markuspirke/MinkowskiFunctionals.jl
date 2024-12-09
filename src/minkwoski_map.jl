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
        ρs = [s.ρ for s in h0_distributions]
        λ = h0_distributions[1].λ
        L = h0_distributions[1].n
        l = floor(Int, L/2)
        Ds = zeros(n - 2l, m - 2l)
        for j in l+1:m-l
            for i in l+1:n-l
                deviation_strengths = zeros(length(ρs))
                for (k, ρ) in enumerate(ρs)
                    bw_map = BWMap(x, ρ)
                    functional = MinkowskiFunctional(bw_map.pixels[i-l:i+l, j-l:j+l])
                    deviation_strengths[k] = deviation_strength(getfield(h0_distributions[k], F), getfield(functional, F))
                end
                Ds[i-l, j-l] = maximum(deviation_strengths)
            end
        end

        return new(λ, ρs, L, Ds)
    end

    # function MinkowskiMap(x::CountsMap, h0_distributions::Vector{MinkowskiDistribution})
    #     m, n = size(x.pixels)
    #     ρs = [s.ρ for s in h0_distributions]
    #     λ = h0_distributions[1].λ
    #     L = h0_distributions[1].n
    #     l = floor(Int, L/2)
    #     Ds = zeros(n - 2l, m - 2l)
    #     Threads.@threads for j in l+1:m-l
    #         deviation_strengths = zeros(length(ρs))
    #         for i in l+1:n-l
    #             for (k, ρ) in enumerate(ρs)
    #                 bw_map = BWMap(x, ρ)
    #                 functional = MinkowskiFunctional(bw_map.pixels[i-l:i+l, j-l:j+l])
    #                 deviation_strengths[k] = deviation_strength(h0_distributions[k], functional)
    #             end
    #             @inbounds Ds[i-l, j-l] = maximum(deviation_strengths)
    #         end
    #     end

    #     return new(λ, ρs, L, Ds)
    # end

    function MinkowskiMap(x::CountsMap, h0_distributions::Vector{MinkowskiDistribution})
        m, n = size(x.pixels)
        ρs = [s.ρ for s in h0_distributions]
        λ = h0_distributions[1].λ
        L = h0_distributions[1].n
        l = floor(Int, L/2)
        Ds = zeros(n - 2l, m - 2l)
        Threads.@threads for j in l+1:m-l
            deviation_strengths = zeros(length(ρs))
            for i in l+1:n-l
                for (k, ρ) in enumerate(ρs)
                    bw_map = BWMap(x, ρ)
                    functional = MinkowskiFunctional(bw_map.pixels[i-l:i+l, j-l:j+l])
                    deviation_strengths[k] = h0_distributions[k].σ[functional]
                end
                @inbounds Ds[i-l, j-l] = maximum(deviation_strengths)
            end
        end

        return new(λ, ρs, L, Ds)
    end

    function MinkowskiMap(x::CountsMap, h0_distributions::Vector{MinkowskiDistribution}, fields)
        m, n = size(x.pixels)
        ρs = [s.ρ for s in h0_distributions]
        λ = h0_distributions[1].λ
        L = h0_distributions[1].n
        l = floor(Int, L/2)
        Ds = zeros(n - 2l, m - 2l)
        Threads.@threads for j in l+1:m-l
            deviation_strengths = zeros(length(ρs))
            for i in l+1:n-l
                for (k, ρ) in enumerate(ρs)
                    bw_map = BWMap(x, ρ)
                    functional = MinkowskiFunctional(bw_map.pixels[i-l:i+l, j-l:j+l])
                    reduced_functional = reduce_functional(functional, fields)
                    deviation_strengths[k] = h0_distributions[k].σ[reduced_functional]
                end
                @inbounds Ds[i-l, j-l] = maximum(deviation_strengths)
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
    mask = probs(h0_distribution) .<= p
    return -log10(sum(probs(h0_distribution)[mask]))
end
