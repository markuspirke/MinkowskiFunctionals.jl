"""
    function MinkowskiMap(x::CountsMap, b::Background, L::Int64)

This calculates a MinkowskiMap for a background gives for example from a FoV background modell
based on only the Area functional. """
function nonhomogeonus_minkowski_map(x::CountsMap, b::Background, L::Int64)
    m, n = size(x.pixels)
    if isodd(L)
        l1 = floor(Int, L/2)
        l2 = floor(Int, L/2)
        αs = ones(n - 2l1, m - 2l2)
        signs = ones(n - 2l1, m - 2l2)
    else
        l1 = floor(Int, (L-1)/2)
        l2 = ceil(Int, (L-1)/2)
        l = l1+l2
        αs = ones(n - l, m - l)
        signs = ones(n - l, m - l)
    end 

    for j in l1+1:m-l2
        for i in l1+1:n-l2
            if b.pixels[i, j] == 0.0
                continue
            end
            local_counts = x[i-l1:i+l2, j-l1:j+l2]
            local_background = b[i-l1:i+l2, j-l1:j+l2]
            ρs = get_thresholds(local_counts)
            l_ρ = length(ρs)
            αs_ρ = ones(l_ρ)
            signs_ρ = zeros(l_ρ)
            for k in 1:length(ρs)
                mink_distribution = AreaDistribution(local_background, ρs[k])
                @inbounds αs_ρ[k] = compatibility(mink_distribution, local_counts)
                @inbounds signs_ρ[k] = get_sign(mink_distribution, local_counts)
            end
            idx = argmin(αs_ρ)
            α = αs_ρ[idx]
            @inbounds αs[i-l1, j-l1] = correct_trials(α, l_ρ)
            @inbounds signs[i-l1, j-l1] = signs_ρ[idx]
        end
    end
    return MinkowskiMap(αs .* signs)
end

"""
    function nonhomogeonus_minkowski_map(x::CountsMap, b::Background, mask::Union{BitMatrix, Matrix{Bool})

This calculates a MinkowskiMap for a background gives for example from a FoV background modell
based on only the Area functional. """
function nonhomogeonus_minkowski_map(x::CountsMap, b::Background, mask::Union{BitMatrix, Matrix{Bool}})
    m, n = size(x.pixels)
    m_mask, n_mask = size(mask)
    l = floor(Int, m_mask/2)
    αs = ones(n - 2l, m - 2l)
    signs = ones(n - 2l, m - 2l)
    for j in l+1:m-l
        for i in l+1:n-l
            if b.pixels[i, j] == 0.0
                continue
            end
            local_counts = x[i-l:i+l, j-l:j+l] .* mask
            local_background = b[i-l:i+l, j-l:j+l] .* mask
            ρs = get_thresholds(local_counts)
            l_ρ = length(ρs)
            αs_ρ = ones(l_ρ)
            signs_ρ = zeros(l_ρ)
            for k in 1:length(ρs)
                mink_distribution = AreaDistribution(local_background, ρs[k])
                @inbounds αs_ρ[k] = compatibility(mink_distribution, local_counts)
                @inbounds signs_ρ[k] = get_sign(mink_distribution, local_counts)
            end
            idx = argmin(αs_ρ)
            α = αs_ρ[idx]
            @inbounds αs[i-l, j-l] = correct_trials(α, l_ρ)
            @inbounds signs[i-l, j-l] = signs_ρ[idx]
        end
    end
    return MinkowskiMap(αs .* signs)
end
