"""
    function lima(Non, Noff, a=1)

This calculates the λ from the Li Ma paper, which is based on the Wilkstheorem.
"""
function λ_lima(Non, Noff, a=1)
    return (a/(a+1) * (Non + Noff)/Non)^Non * (1/(a+1) * (Non + Noff)/Noff)^Noff
end

function significance_lima_exact(Non, λ)
    x = Non * log(Non/λ) + λ - Non
    if x < 0.0
        return 0.0
    else
        return sqrt(2) * sqrt(x)
    end
end


"""
    function significance_lima(Non, Noff, a=1)

This function transforms the λ from Wilks theorem and return a random variable which
is distributed according to a standard normal distribution.
That means we can directely read of the significance from this value.
"""
function significance_lima(Non, Noff, a=1)
    return sqrt(-2*log(λ_lima(Non, Noff, a)))
end

"""
    function lima_map(img::CountsMap, background::Background, L::Int64)

This takes a counts map an returns a significance map based on the
formula from the lima paper. It uses an extended region to calculate
the signifiance of the center pixels by summing up all counts in the region,
these are the ON counts and as OFF counts the correct λ is used which is
multiplied by the number of pixels L^2, where L is the window size.
"""
function lima_map(img::CountsMap, background::Background, L::Int64)
    m, n = size(img)
    l = floor(Int, L/2)
    r = floor(Int, L/2)
    σs = zeros(n - 2l, m - 2l)

    for j in l+1:m-l
        for i in l+1:n-l
            window = img[i-l:i+l, j-l:j+l]# .* round_mask
            b_window = background[i-l:i+l, j-l:j+l]# .* round_mask
            σ = significance_lima_exact.(sum(window), sum(b_window))
            s = sum(window) > sum(b_window) ? 1.0 : -1.0
            σs[i-l,j-l] = σ * s
        end
    end
    σs
end

"""
    function lima_map(img::CountsMap, background::Matrix{Float64}, L::Int64)

This takes a counts map an returns a significance map based on the
formula from the lima paper. It uses an extended region to calculate
the signifiance of the center pixels by summing up all counts in the region,
these are the ON counts and as OFF counts the correct λ is used which is
multiplied by the number of pixels L^2, where L is the window size.
"""
function lima_map(img::CountsMap, background::Float64, L::Int64)
    m, n = size(img)
    l = floor(Int, L/2)
    r = floor(Int, L/2)
    σs = zeros(n - 2l, m - 2l)
    for j in l+1:m-l
        for i in l+1:n-l
            window = img[i-l:i+l, j-l:j+l]# .* round_mask
            σ = significance_lima_exact.(sum(window), L^2*background)
            s = sum(window) > L^2*background ? 1.0 : -1.0
            σs[i-l,j-l] = σ * s
        end
    end
    σs
end
"""
    function lima_map(img::CountsMap, background::Matrix{Float64}, mask::Matrix{Bool})

This takes a counts map an returns a significance map based on the
formula from the lima paper. It uses an extended region to calculate
the signifiance of the center pixels by summing up all counts in the region,
these are the ON counts and as OFF counts the correct λ is used which is
multiplied by the number of pixels L^2, where L is the window size.
The difference to the above method is that here we do not use a square
window but a round one.
"""
function lima_map(img::CountsMap, background::Background, mask::Union{BitMatrix, Matrix{Bool}})
    m, n = size(img)
    m_mask, n_mask = size(mask)
    l = floor(Int, m_mask/2)
    σs = zeros(n - 2l, m - 2l)

    for j in l+1:m-l
        for i in l+1:n-l
            window = img[i-l:i+l, j-l:j+l] .* mask
            b_window = background[i-l:i+l, j-l:j+l] .* mask
            σ = significance_lima_exact.(sum(window), sum(b_window))
            s = sum(window) > sum(b_window) ? 1.0 : -1.0
            σs[i-l,j-l] = σ * s
        end
    end
    σs
end
