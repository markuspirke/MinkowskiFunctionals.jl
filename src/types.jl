struct CountsMap
    pixels
end

function CountMap(d::Poisson, n)
    return CountMap(rand(d, n, n))
end

function CountMap(d::FullNormal, n, N)
    samples = rand(d, N)
    counts_map = fit(Histogram, (samples[1, :], samples[2,:]), (range(-1, 1, n+1), range(-1,1,n+1)))

    return CountMap(counts_map.weights)
end

struct BWMap
    pixels
end

function BWMap(x::CountsMap, ρ)
    return BWMap(x.counts .> ρ)
end

function rand(T::Type{BWMap}, n, m)
    return T(reshape(rand(Bool, n*m), n, m))
end

struct MinkowskiFunctional
    A
    P
    χ
end

struct BasesState
    pixels::SMatrix{2, 2}
    functional::MinkowskiFunctional
end
