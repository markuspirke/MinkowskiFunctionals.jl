"""
This generates a Minkowski map from a given CountsMap.
This also needs distributions of Minkowski functionals
and one needs to decide which functional is used.
"""
struct MinkowskiMap
    pixels::Matrix{Float64}
end

function p2σ(x::MinkowskiMap)
    p2σ.(abs.(x.pixels)) .* sign.(x.pixels)
end

Base.size(x::MinkowskiMap) = size(x.pixels)
Base.getindex(x::MinkowskiMap, i, j) = x.pixels[i, j]

get_sign(d::AreaDistribution, x::CountsMap) = d.p.p*d.p.n > sum(BWMap(x, d.ρ).pixels) ? -1.0 : 1.0
get_sign(d::AreaDistribution, x::Matrix{Int64}) = d.p.p*d.p.n > sum(BWMap(x, d.ρ).pixels) ? -1.0 : 1.0


"""
    function correct_trials(x::Float64, n::Int64)

If significances are calculated multiply times, by looking differently on the data, trial factors are needed.
This function calculates trial factors exactly for p values greater than 1e-16 and approximates them with
a Taylor expansion for smaller values.
"""
function correct_trials(x::Float64, n::Int64)
    if x > 1e-16
        return 1 - (1 - x)^n
    else
        return n*x
    end
end

"""
    function correction!(x, b, target)

This accounts for changing background within the local window and either adds or subtracts counts.
If the background at a particular pixel is smaller than the background of the center pixel
the counts are filled up by sampling from a Poisson distribution. If otherwise, so less counts at the particular
pixel counts are removed by sampling from a binomial distribution.
"""
function correction!(x, b, target)
    m, n = size(x)
    for j in 1:n
        for i in 1:m
            if target - b[i, j] != 0.0
                sign_corr = sign(target - b[i,j])
                correction = sign_corr > 0.0 ? rand(Poisson(abs(target - b[i,j]))) : -1.0 * rand(Binomial(x[i, j], target/b[i, j]))
                x[i, j] += correction
                if x[i, j] < 0
                    x[i, j] = 0
                end
            end
        end
    end
end

"""
    function get_tresholds(local_counts, local_background)

Calculates which tresholds are used.
"""
function get_tresholds(local_counts, local_background)
    ρ_min = minimum(local_counts) > 0 ? minimum(local_counts) : 1
    ρ_max = maximum(local_counts) > 0 ? maximum(local_counts) : 1
    return ρ_min:ρ_max
end

"""
    function MinkowskiMap(x::CountsMap, b::Background, L::Int64)

This calculates a MinkowskiMap for a background gives for example from a FoV background modell
based on only the Area functional.
"""
function MinkowskiMap(x::CountsMap, b::Background, L::Int64)
    m, n = size(x.pixels)
    l = floor(Int, L/2)
    αs = ones(n - 2l, m - 2l)
    signs = ones(n - 2l, m - 2l)
    for j in l+1:m-l
        for i in l+1:n-l
            if b.pixels[i, j] == 0.0
                continue
            end
            local_counts = x[i-l:i+l, j-l:j+l]
            local_background = b[i-l:i+l, j-l:j+l]
            # if i == 218 && j == 51
            #     @show local_counts
            #     @show local_background
            # end
            correction!(local_counts, local_background, b[i, j])
            ρs = get_tresholds(local_counts, local_background)
            l_ρ = length(ρs)
            αs_ρ = ones(l_ρ)
            signs_ρ = zeros(l_ρ)
            Threads.@threads for k in 1:length(ρs)
                mink_distribution = AreaDistribution(L^2, b.pixels[i, j], ρs[k])
                αs_ρ[k] = compatibility(mink_distribution, x[i-l:i+l, j-l:j+l])
                # if αs_ρ[k] < 5e-7
                #     @show i, j
                #     @show local_counts
                #     @show local_background
                #     return
                # end
                signs_ρ[k] = get_sign(mink_distribution, local_counts)
            end
            idx = argmin(αs_ρ)
            α = αs_ρ[idx]
            @inbounds αs[i-l, j-l] = correct_trials(α, l_ρ)
            @inbounds signs[i-l, j-l] = signs_ρ[idx]
        end
    end
    return MinkowskiMap(αs .* signs)
end

"""
    function MinkowskiMap(x::CountsMap, b::Background, mask::Union{BitMatrix, Matrix{Bool}})

This calculates a MinkowskiMap for a background gives for example from a FoV background modell
based on only the Area functional. A different kernel can be given as a mask.
"""
function MinkowskiMap(x::CountsMap, b::Background, mask::Union{BitMatrix, Matrix{Bool}})
    m, n = size(x.pixels)
    m_mask, n_mask = size(mask)
    N = sum(mask)
    l = floor(Int, m_mask/2)
    αs = ones(n - 2l, m - 2l)
    signs = ones(n - 2l, m - 2l)
    for j in l+1:m-l
        for i in l+1:n-l
            if b.pixels[i, j] == 0.0
                continue
            end
            local_counts = x[i-l:i+l, j-l:j+l]
            local_background = b[i-l:i+l, j-l:j+l]
            correction!(local_counts, local_background, b[i, j])
            ρs = get_tresholds(local_counts, local_background)
            l_ρ = length(ρs)
            αs_ρ = ones(l_ρ)
            signs_ρ = zeros(l_ρ)
            Threads.@threads for k in 1:length(ρs)
                mink_distribution = AreaDistribution(N, b.pixels[i, j], ρs[k])
                αs_ρ[k] = compatibility(mink_distribution, x[i-l:i+l, j-l:j+l] .* mask)
                signs_ρ[k] = get_sign(mink_distribution, local_counts)
            end
            idx = argmin(αs_ρ)
            α = αs_ρ[idx]
            @inbounds αs[i-l, j-l] = correct_trials(α, l_ρ)
            @inbounds signs[i-l, j-l] = signs_ρ[idx]
        end
    end
    return MinkowskiMap(αs .* signs)
end

"""
    function MinkowskiMap(x::CountsMap, b::Float64, L::Int64)

This calculates a MinkowskiMap for a fixed and constant background based on only the Area functional.
"""
function MinkowskiMap(x::CountsMap, b::Float64, L::Int64)

    m, n = size(x.pixels)
    ρs = b:maximum(x.pixels)
    l_ρ = length(ρs)
    l = floor(Int, L/2)
    αs = zeros(n - 2l, m - 2l)
    signs = zeros(n - 2l, m - 2l)
    mink_ds = [AreaDistribution(L^2, b, ρ) for ρ in ρs]
    for j in l+1:m-l
        for i in l+1:n-l
            αs_ρ = zeros(l_ρ)
            signs_ρ = zeros(l_ρ)
            for (k, ρ) in enumerate(ρs)
                αs_ρ[k] = compatibility(mink_ds[k], x[i-l:i+l, j-l:j+l])
                signs_ρ[k] = get_sign(mink_ds[k], x[i-l:i+l, j-l:j+l])
            end
            idx = argmin(αs_ρ)
            α = αs_ρ[idx]
            @inbounds αs[i-l, j-l] = 1 - (1 - α)^l_ρ
            @inbounds signs[i-l, j-l] = signs_ρ[idx]
        end
    end

    return MinkowskiMap(αs .* signs)
end

"""
    function MinkowskiMap(x::CountsMap, b::Background, L::Int64)

This calculates a MinkowskiMap for a background gives for example from a FoV background modell
based on the information by all functionals.
"""
function MinkowskiMap(x::CountsMap, b::Background, Ω::DensityOfStates)
    m, n = size(x.pixels)
    L = Ω.n
    l = floor(Int, L/2)
    αs = ones(n - 2l, m - 2l)
    signs = ones(n - 2l, m - 2l)
    for j in l+1:m-l
        for i in l+1:n-l
            if b.pixels[i, j] == 0.0
                continue
            end
            local_counts = x[i-l:i+l, j-l:j+l]
            local_background = b[i-l:i+l, j-l:j+l]
            correction!(local_counts, local_background, b[i, j])
            ρs = get_tresholds(local_counts, local_background)
            l_ρ = length(ρs)
            αs_ρ = ones(l_ρ)
            signs_ρ = zeros(l_ρ)
            # for (k, ρ) in enumerate(ρs)
            # for k in 1:length(ρs)
            Threads.@threads for k in 1:length(ρs)
                mink_distribution = MinkowskiDistribution(Ω, b.pixels[i, j], ρs[k])
                αs_ρ[k] = compatibility(mink_distribution, x[i-l:i+l, j-l:j+l])
                signs_ρ[k] = get_sign(mink_distribution, local_counts)
            end
            idx = argmin(αs_ρ)
            α = αs_ρ[idx]
            @inbounds αs[i-l, j-l] = correct_trials(α, l_ρ)
            @inbounds signs[i-l, j-l] = signs_ρ[idx]
        end
    end
    return MinkowskiMap(αs .* signs)
end


"""
    function MinkowskiMap(x::CountsMap, b::Float64, L::Int64)

This calculates a MinkowskiMap for a fixed and constant background
based on the information by all functionals.
"""
function MinkowskiMap(x::CountsMap, b::Float64, Ω::DensityOfStates)
    m, n = size(x.pixels)

    ρs = b:maximum(x.pixels)
    l_ρ = length(ρs)
    L = Ω.n
    l = floor(Int, L/2)
    αs = zeros(n - 2l, m - 2l)
    signs = zeros(n - 2l, m - 2l)
    mink_ds = [MinkowskiDistribution(Ω, b, ρ) for ρ in ρs]
    for j in l+1:m-l
        for i in l+1:n-l
            αs_ρ = zeros(l_ρ)
            signs_ρ = zeros(l_ρ)
            Threads.@threads for k in 1:length(ρs)
                αs_ρ[k] = compatibility(mink_ds[k], x[i-l:i+l, j-l:j+l])
                signs_ρ[k] = get_sign(mink_ds[k], x[i-l:i+l, j-l:j+l])
            end
            idx = argmin(αs_ρ)
            α = αs_ρ[idx]
            @inbounds αs[i-l, j-l] = 1 - (1 - α)^l_ρ
            @inbounds signs[i-l, j-l] = signs_ρ[idx]
        end
    end

    return MinkowskiMap(αs .* signs)
end

get_sign(d::MinkowskiDistribution, x::CountsMap) = d.p_black*d.n^2 > sum(BWMap(x, d.ρ).pixels) ? -1.0 : 1.0
get_sign(d::MinkowskiDistribution, x::Matrix{Int64}) = d.p_black*d.n^2 > sum(BWMap(x, d.ρ).pixels) ? -1.0 : 1.0
