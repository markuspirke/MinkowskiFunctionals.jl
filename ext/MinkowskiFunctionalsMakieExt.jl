module MinkowskiFunctionalsMakieExt

using MinkowskiFunctionals, CairoMakie
isdefined(Base, :get_extension) ? (using CairoMakie) : (using ..CairoMakie)

function CairoMakie.heatmap(bw_map::BWMap; kwargs...)
    fig = Figure()
    ax = Axis(fig[1, 1], aspect = DataAspect())
    hm = heatmap!(ax, bw_map.pixels .== 0, colormap=:grays; kwargs...)
    hidedecorations!(ax)

    fig, ax, hm
end

function CairoMakie.heatmap(counts_map::CountsMap; kwargs...)
    fig = Figure()
    ax = Axis(fig[1, 1], aspect = DataAspect())
    hm = heatmap!(ax, counts_map.pixels, colormap=:afmhot; kwargs...)
    hidedecorations!(ax)

    fig, ax, hm
end

function CairoMakie.heatmap(mink_map::MinkowskiMap; kwargs...)
    fig = Figure()
    ax = Axis(fig[1, 1], aspect = DataAspect())
    hm = heatmap!(ax, p2σ(mink_map), colormap=:vik; kwargs...)
    Colorbar(fig[1, 2], hm, label="σ")
    hidedecorations!(ax)

    fig, ax, hm
end

function CairoMakie.hist(mink_map::MinkowskiMap; kwargs...)
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel="σ")
    h = hist!(ax, p2σ(mink_map)[1:end]; kwargs...)

    fig, ax, h
end

end
