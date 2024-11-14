# Introduction

This is a package which calculates Minkowski functional in two dimensions. The package is written as a tool for $\gamma$-ray astronomy, but is in principal generic.

## Calculate Minkowski functionals

First let us take a look how to compute the Minkowski functionals for some random black and white image.
```julia
using MinkowskiFunctionals
img = rand(Bool, 10, 10)
result = MinkowskiFunctional(img)
```

## Exact distributions of Poisson random fields

With a brute force approach we can calculate the exact distributions of Minkowski functionals for Poisson random fields up to size 6 x 6 pixels. This can be done like this
```julia
using MinkowskiFunctionals

n = 3
λ = 10
ρ = 10

exact_distributions = PoissonMinkowskiDistributions(n, λ, ρ)
```
We can access the individual distribution for the area, perimeter and Euler characterisic like this
```julia
exact_distributions.A
exact_distributions.P
exact_distributions.χ
```

## Sampling distributions for Poisson random fields
Even for 5 x 5 pixels generating exact distributions is computationally heavy. Instead we can also sample distributions like this

```julia
N = 10000
n = 3
λ = 10
ρ = 10
sampled_distributions = SampledPoissonMinkowskiDistributions(N, n, λ, ρ)

sampled_distributions.A
sampled_distributions.P
sampled_distributions.χ
```
