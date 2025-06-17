### A Pluto.jl notebook ###
# v0.17.0

using Markdown
using InteractiveUtils

# ╔═╡ a98fae2c-9c3a-41c6-96f3-93d147f79e7b
begin
    using Pkg
    Pkg.activate(".testenv")
    Pkg.add("Revise")
    using Revise
    Pkg.add(["PlutoUI", "ExtendableGrids", "PyPlot", "Plots", "GLMakie"])
    Pkg.develop(["GridVisualize", "PlutoVista"])
end

# ╔═╡ 9701cbe0-d048-11eb-151b-67dda7b72b71
begin
    using PlutoVista
    using PyPlot
    using GLMakie
    using Plots
    using GridVisualize
    using ExtendableGrids
    using PlutoUI
end

# ╔═╡ 4ab54a93-d759-4d85-9d74-47a220afaa2c
PyPlot.svg(true)

# ╔═╡ 68e2c958-b417-4ba1-9577-697304fe140a
TableOfContents()

# ╔═╡ b35d982d-1fa9-413d-b008-892b4f241097
md"""
# Scalarplot2D test
"""

# ╔═╡ d5258595-60e4-406f-a71e-69111cdad8b9
function testplot_rect(; Plotter = nothing, kwargs...)
    X = 0:0.05:5
    grid = simplexgrid(X, X)
    f = map((x, y) -> sinpi(x) * cos(x * y), grid)
    return scalarplot(grid, f; Plotter = Plotter, kwargs...)
end

# ╔═╡ 0998a9a7-d57a-476e-aacd-bee9396e9b8f
testplot_rect(; Plotter = PyPlot, colorbarticks = 15, resolution = (500, 300))

# ╔═╡ 55cc2e7b-bb1e-437d-9377-75142bb003ec
testplot_rect(;
    Plotter = PlutoVista, colorbarticks = 11, resolution = (400, 300),
    backend = :plotly, colormap = :hot, levels = 0
)

# ╔═╡ a8cdce79-bd03-40d7-9b4d-5336b4bb1056
testplot_rect(;
    Plotter = PlutoVista, colorbarticks = 15, resolution = (400, 300),
    backend = :vtk, colormap = :hot
)

# ╔═╡ 94a905f9-83dc-47a5-8a4e-fa275ffc3408
testplot_rect(; Plotter = GLMakie, colorbarticks = 15, resolution = (300, 300))

# ╔═╡ 00ae8154-7753-4374-a0e9-ead8ed572b24
testplot_rect(; Plotter = Plots, levels = 5, colorbarticks = 20, resolution = (400, 300))

# ╔═╡ ba5111b8-0dca-42d2-970f-1e88f5392324
html"""<hr>"""

# ╔═╡ 92bfccf7-abf1-47d5-8d8b-9ae9003ad1ac
md"""
## Appendix
"""

# ╔═╡ Cell order:
# ╠═9701cbe0-d048-11eb-151b-67dda7b72b71
# ╠═4ab54a93-d759-4d85-9d74-47a220afaa2c
# ╠═68e2c958-b417-4ba1-9577-697304fe140a
# ╟─b35d982d-1fa9-413d-b008-892b4f241097
# ╠═d5258595-60e4-406f-a71e-69111cdad8b9
# ╠═0998a9a7-d57a-476e-aacd-bee9396e9b8f
# ╠═55cc2e7b-bb1e-437d-9377-75142bb003ec
# ╠═a8cdce79-bd03-40d7-9b4d-5336b4bb1056
# ╠═94a905f9-83dc-47a5-8a4e-fa275ffc3408
# ╠═00ae8154-7753-4374-a0e9-ead8ed572b24
# ╟─ba5111b8-0dca-42d2-970f-1e88f5392324
# ╟─92bfccf7-abf1-47d5-8d8b-9ae9003ad1ac
# ╠═a98fae2c-9c3a-41c6-96f3-93d147f79e7b
