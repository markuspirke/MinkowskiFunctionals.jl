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

function CairoMakie.plot(bw_map, f::MinkowskiFunctional)
    n, m = size(bw_map)
    img = zeros(n+2, n+2)
	img[2:end-1, 2:end-1] = bw_map
    bw_map = img .== 1.0

    fig = Figure()
    ax1 = Axis(fig[1, 1],
            xticks=1:n+2, yticks=1:m+2,
            xgridcolor = :gray,
            ygridcolor = :gray,
            xgridwidth = 2,
            ygridwidth = 2,
            title = "Black and White Image",
            aspect=DataAspect())

    ax2 = Axis(fig[1, 2],
            xticks=1:n+2, yticks=1:m+2,
            xgridcolor = :gray,
            ygridcolor = :gray,
            xgridwidth = 2,
            ygridwidth = 2,
            title = "Area A=$(f.A)",
            aspect=DataAspect())

    ax3 = Axis(fig[2, 1],
            xticks=1:n+2, yticks=1:m+2,
            xgridcolor = :gray,
            ygridcolor = :gray,
            xgridwidth = 2,
            ygridwidth = 2,
            title = "Perimeter P=$(f.P)",
            aspect=DataAspect())

    ax4 = Axis(fig[2, 2],
            xticks=1:n+2, yticks=1:m+2,
            xgridcolor = :gray,
            ygridcolor = :gray,
            xgridwidth = 2,
            ygridwidth = 2,
            title = "Euler Characteristic Χ:$(f.χ)",
            aspect=DataAspect())

    hm = heatmap!(ax1, 0:n+2, 0:m+2, bw_map .== 0, colormap=:grays)
    translate!(hm, 0, 0, -100)
    hidedecorations!(ax1, grid = false)

    hm = heatmap!(ax2, 0:n+2, 0:m+2, bw_map, colormap=[:white, :lightcoral])
    translate!(hm, 0, 0, -100)
    hidedecorations!(ax2, grid = false)

    edges = MinkowskiFunctionals.find_perimeter_edges(bw_map)
    hm = heatmap!(ax3, 0:n+2, 0:n+2, bw_map .== 0, colormap=:grays)
	for e in edges
		lines!(ax3, [e.a, e.b],
			   color=:cornflowerblue, linewidth=2, linecap = :square)
	end
	translate!(hm, 0, 0, -100)
	hidedecorations!(ax3, grid = false)

    intermed_bw_map = MinkowskiFunctionals.label_outer_background(bw_map)
	outer_edges = MinkowskiFunctionals.find_perimeter_edges(intermed_bw_map)
	inner_edges = MinkowskiFunctionals.Edge[]
	for e in edges
		e in outer_edges ? continue : push!(inner_edges, e)
	end
    hm = heatmap!(ax4, 0:n+2, 0:n+2, bw_map .== 0, colormap=:grays)
	for e in outer_edges
		lines!(ax4, [e.a, e.b],
			   color=:plum, linewidth=2, linecap = :square)
	end
	for e in inner_edges
		lines!(ax4, [e.a, e.b],
			   color=:gold, linewidth=2, linecap = :square)
	end
	translate!(hm, 0, 0, -100)
	hidedecorations!(ax4, grid = false)

    return fig
end

function CairoMakie.plot(bw_map::BWMap, f::MinkowskiFunctional)
    CairoMakie.plot(bw_map.pixels, f)
end
end
