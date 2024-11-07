using Documenter, Pkg
using Pkg
Pkg.develop(path="..")
using MinkowskiFunctionals

makedocs(
    sitename="MinkowskiFunctionals.jl",
    authors="Markus Pirke",
)
