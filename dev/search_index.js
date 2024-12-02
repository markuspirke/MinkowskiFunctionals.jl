var documenterSearchIndex = {"docs":
[{"location":"#Introduction","page":"Introduction","title":"Introduction","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"This is a package which calculates Minkowski functional in two dimensions. The package is written as a tool for gamma-ray astronomy, but is in principal generic.","category":"page"},{"location":"#Calculate-Minkowski-functionals","page":"Introduction","title":"Calculate Minkowski functionals","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"First let us take a look how to compute the Minkowski functionals for some random black and white image.","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"using MinkowskiFunctionals\nimg = rand(Bool, 10, 10)\nresult = MinkowskiFunctional(img)","category":"page"},{"location":"#Exact-distributions-of-Poisson-random-fields","page":"Introduction","title":"Exact distributions of Poisson random fields","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"With a brute force approach we can calculate the exact distributions of Minkowski functionals for Poisson random fields up to size 6 x 6 pixels. This can be done like this","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"using MinkowskiFunctionals\n\nn = 3\nλ = 10\nρ = 10\n\nexact_distributions = PoissonMinkowskiDistributions(n, λ, ρ)","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"We can access the individual distribution for the area, perimeter and Euler characterisic like this","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"exact_distributions.A\nexact_distributions.P\nexact_distributions.χ","category":"page"},{"location":"#Sampling-distributions-for-Poisson-random-fields","page":"Introduction","title":"Sampling distributions for Poisson random fields","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"Even for 5 x 5 pixels generating exact distributions is computationally heavy. Instead we can also sample distributions like this","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"N = 10000\nn = 3\nλ = 10\nρ = 10\nsampled_distributions = SampledPoissonMinkowskiDistributions(N, n, λ, ρ)\n\nsampled_distributions.A\nsampled_distributions.P\nsampled_distributions.χ","category":"page"}]
}