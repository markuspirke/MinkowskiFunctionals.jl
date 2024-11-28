"""
    function deviation2σ(D)

Converts deviation strenght to a significance as given in σ of a unit Gaussian.
"""
function deviation2σ(D)
	p = 1 - 10^(-D)
	return quantile(Normal(0, 1), 1/2+p/2) # till p=0.5 we have 0, then go for half of p again
end
