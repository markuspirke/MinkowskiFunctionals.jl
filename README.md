# MinkowskiFunctionals

[![Build Status](https://github.com/markuspirke/MinkowskiFunctionals.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/markuspirke/MinkowskiFunctionals.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://markuspirke.github.io/MinkowskiFunctionals.jl/dev/)

**MinkowskiFunctionals.jl** is a Julia package for calculating 2D Minkowski functionals.

Minkowski functionals offer a powerful method for characterizing the shape, size, and connectivity of structures in binary images, without relying on overly detailed descriptions and assumptions.

This package allows users to compute the functionals area, perimeter, and Euler characteristic. It can generate the exact distributions with an "brute force approach", by essentially going through all possible black and white images (not advised to try out on images larger than 6 x 6 pixels). 

If the exact distribution is not necessary, the package provides functionality to approximately sample these distributions.

## Installation Guide

1. Download the latest Julia version ([latest release](https://julialang.org/downloads/) or [via juliaup](https://github.com/JuliaLang/juliaup))
2. Open a terminal and hit **julia**. This will open a **Julia REPL**.
3. Inside the REPL hit ]. This will open up the Julia package manager. Then type **add MinkowskiFunctionals**.
```julia
julia>]
(v1.8) pkg> add MinkwoskiFunctionals
```

## Basic Usage

```julia
using MinkowskiFunctionals
img = rand(Bool, 10, 10)
result = MinkowskiFunctional(img)
```
