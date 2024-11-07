# MinkowskiFunctionals

[![Build Status](https://github.com/markuspirke/MinkowskiFunctionals.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/markuspirke/MinkowskiFunctionals.jl/actions/workflows/CI.yml?query=branch%3Amain)

**MinkowskiFunctionals.jl** is a Julia package for calculating 2D Minkowski functionals on black-and-white images.
Minkowski functionals offer a powerful method for characterizing the shape, size, and connectivity of structures in binary images, without relying on overly detailed descriptions and assumptions.
This package allows users to compute the functionals area, perimeter, and Euler characteristic. It can generate the exact distributions with an "brute force approach", by essentially going through all possible black and white images (not advised to try out on images larger than 7 x 7 pixels). 
If the exact distribution is not necessary, the package provides functionality to approximately sample these distributions.

## Basic Usage

``` julia
using MinkowskiFunctionals
img = rand(Bool, 10, 10)
result = MinkowskiFunctional(img)
```
