# Introduction

This is a package which calculates Minkowski functional in two dimensions. The package is written as a tool for $\gamma$-ray astronomy, but is in principal generic.

## Usage
First let us take a look how to compute the Minkowski functionals for some random black and white image.
``` julia
using MinkowskiFunctionals
img = rand(Bool, 10, 10)
result = MinkowskiFunctional(img)
```
