"""
    function lima(Non, Noff, a=1)

This calculates the λ from the Li Ma paper, which is based on the Wilkstheorem.
"""
function λ_lima(Non, Noff, a=1)
    return (a/(a+1) * (Non + Noff)/Non)^Non * (1/(a+1) * (Non + Noff)/Noff)^Noff
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
    function lima_map(img, L, λ)

This takes a counts map an returns a significance map based on the
formula from the lima paper. It uses an extended region to calculate
the signifiance of the center pixels by summing up all counts in the region,
these are the ON counts and as OFF counts the correct λ is used which is
multiplied by the number of pixels L^2, where L is the window size.
"""
function lima_map(img::CountsMap, L, λ)
    m, n = size(img)
    l = floor(Int, L/2)
    σs = zeros(n - 2l, m - 2l)
    for j in l+1:m-l
        for i in l+1:n-l
            σ = significance_lima.(sum(img[i-l:i+l, j-l:j+l]), L^2*λ)
            σs[i-l,j-l] = σ
        end
    end
    σs
end

"""
    function lima_map_roundkernel(img, L, λ)

This takes a counts map an returns a significance map based on the
formula from the lima paper. It uses an extended region to calculate
the signifiance of the center pixels by summing up all counts in the region,
these are the ON counts and as OFF counts the correct λ is used which is
multiplied by the number of pixels L^2, where L is the window size.
The difference to the above method is that here we do not use a square
window but a round one.
"""
function lima_map_roundkernel(img::CountsMap, L, λ)
    m, n = size(img)
    l = floor(Int, L/2)
    r = floor(Int, L/2)
    σs = zeros(n - 2l, m - 2l)
    round_mask = [sqrt((i - r - 1)^2 + (j - r - 1)^2) <= r for i in 1:L, j in 1:L]

    for j in l+1:m-l
        for i in l+1:n-l
            window = img[i-l:i+l, j-l:j+l]
            round_window = window .* round_mask
            σ = significance_lima.(sum(round_window), sum(round_mask)*λ)
            σs[i-l,j-l] = σ
        end
    end
    σs
end
