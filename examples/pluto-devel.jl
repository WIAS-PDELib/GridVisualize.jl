### A Pluto.jl notebook ###
# v0.16.1

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    return quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 6df3beed-24a7-4b26-a315-0520f4863190
develop = true

# ╔═╡ 9701cbe0-d048-11eb-151b-67dda7b72b71
begin
    using Pkg
    Pkg.activate(mktempdir())
    Pkg.add("Revise")
    using Revise
    Pkg.add(["PyPlot", "PlutoUI", "Plots"])
    if develop
        Pkg.add("GLMakie")
        Pkg.develop(["PlutoVista", "GridVisualize", "ExtendableGrids"])
    else
        Pkg.add(["PlutoVista", "GridVisualize", "ExtendableGrids"])
    end
    using PlutoVista, GridVisualize, ExtendableGrids, PlutoUI
    import Plots, PyPlot
    if develop
        import GLMakie
    else
        GLMakie = nothing
    end
end

# ╔═╡ 62c6c85c-e374-4315-b737-b7a0005d2985
Plots.gr()

# ╔═╡ 68e2c958-b417-4ba1-9577-697304fe140a
TableOfContents()

# ╔═╡ b35d982d-1fa9-413d-b008-892b4f241097
md"""
# Test notebook for PlutoVista backend of GridVisualize
"""

# ╔═╡ 00b04f6b-34a6-4e30-9864-d273305281d4
md"""
## 1D scalar plot
"""

# ╔═╡ a20d74c9-16da-408a-b247-0c17321888f9
function testplot1(Plotter)
    grid = simplexgrid(0:0.01:10)
    return scalarplot(
        grid, map(sin, grid); Plotter = Plotter, resolution = (600, 200),
        markershape = :star5, markevery = 20, xlabel = "x", ylabel = "z",
        legend = :rt, label = "sin"
    )
end

# ╔═╡ cc17187f-404c-4c31-8625-fa067eea7273
testplot1(PyPlot)

# ╔═╡ ad1c4dd6-7c6b-4433-a1ac-b2f817ba5d81
testplot1(Plots)

# ╔═╡ 33482af8-3542-4723-ae43-770a789b69b3
testplot1(PlutoVista)

# ╔═╡ d21982eb-aa7e-424c-932e-f0b1460f19b9
testplot1(GLMakie)

# ╔═╡ c4eeb06f-932e-4acc-8e5b-f2a7f9242a42
function testplot2(Plotter; t = 0)
    p = GridVisualizer(;
        Plotter = Plotter, resolution = (500, 200), legend = :rt,
        xlabel = "x"
    )
    grid = simplexgrid(0:0.01:10)
    scalarplot!(
        p, grid, map(x -> sin(x - t), grid); Plotter = Plotter, color = :red,
        label = "sin(x-$(t))", linestyle = :dash
    )
    scalarplot!(
        p, grid, map(cos, grid); Plotter = Plotter, color = :green, clear = false,
        label = "cos", linestyle = :dashdot, linewidth = 3
    )
    return reveal(p)
end

# ╔═╡ 29ca4775-6ba5-474c-bd2c-8f770b09addd
testplot2(PyPlot)

# ╔═╡ c0ab77e8-01ea-436d-85f2-34e253944f11
testplot2(Plots)

# ╔═╡ 4f73fb99-fab9-4ca4-addb-7d12d16ceba0
testplot2(GLMakie)

# ╔═╡ 84192945-d4b6-4949-8f06-d94e04a7a56d
testplot2(PlutoVista)

# ╔═╡ 7fbaf93f-3cfb-47d0-8252-487e60ba3e54
md"""
### Changing data by re-creating a plot
"""

# ╔═╡ 63fe3259-7d79-40ec-98be-e0592e40ee6b
@bind t2 PlutoUI.Slider(0:0.1:5, show_value = true)

# ╔═╡ 4de6b5c9-4d2d-4bcb-bc88-c6f50a23f9b6
testplot2(PlutoVista; t = t2)

# ╔═╡ 2061e7fd-c740-4d4b-af5b-7a3a9444aafd
md"""
### Changing data by updating the plot

For this pattern, we observe a notable difference to other backends:  the plot with PlutoVista appears above the cell where the GridVisualizer is created instead of the cell where data are plotted. The reason is that for updating data, we need to have a visualization context which stays the same.
"""

# ╔═╡ f84beb4f-4136-4e5a-ba43-279b703fc75f
begin
    X2 = 0:0.001:10
    grid2 = simplexgrid(collect(X2))
    f2(t) = map((x) -> sin(x^2 - t), grid2)
end

# ╔═╡ 29fa4467-65ee-4dad-a660-5197864ddbdc
md"""
t4: $(@bind t4 PlutoUI.Slider(-10:0.1:10, default=0, show_value=true))
"""

# ╔═╡ c1278fb2-3e75-445f-893a-b8b8a7e931d3
p = GridVisualizer(; Plotter = PlutoVista, resolution = (600, 200), dim = 1, legend = :lt);
p;

# ╔═╡ 661531f7-f740-4dd4-9a59-89ddff06ba5c
scalarplot!(p, X2, f2(t4); show = true, clear = true, color = :red, label = "t=$(t4)")

# ╔═╡ ed9b80e5-9678-4ba6-bb36-c2e0674ed9ba
md"""
## 1D grid plot 
"""

# ╔═╡ 9ce4f63d-cd96-48d7-a637-07cb84fa88ab
function testgridplot(Plotter)
    grid = simplexgrid(0:1:10)
    cellmask!(grid, [0.0], [5], 2)
    bfacemask!(grid, [5], [5], 3)
    return gridplot(grid; Plotter = Plotter, resolution = (600, 200), legend = :rt)
end

# ╔═╡ 77eeefc7-e416-426b-8f87-1bc8439dae6d
testgridplot(PyPlot)

# ╔═╡ d503ee1e-1e1f-4235-b286-dc3137a2c96a
testgridplot(PlutoVista)

# ╔═╡ 395d05f4-46a9-4e85-804c-35c79acc96e3
testgridplot(GLMakie)

# ╔═╡ ae1fe1ab-4a0e-4c80-bd6f-912201fb4bb4
md"""
## 2D Scalar plot
"""

# ╔═╡ d5258595-60e4-406f-a71e-69111cdad8b9
function testplot3(Plotter)
    X = 0:0.1:10
    grid = simplexgrid(X, X)
    f = map((x, y) -> sin(x) * atan(y), grid)
    return scalarplot(
        grid, f; Plotter = Plotter,
        resolution = (300, 300), limits = (-π / 2, π / 2)
    )
end

# ╔═╡ c98a90bf-1a3e-4681-a3b0-663c6844df6b
testplot3(PyPlot)

# ╔═╡ a0c3067b-3aa5-493e-b132-89746483b5ce
testplot3(Plots)

# ╔═╡ 0998a9a7-d57a-476e-aacd-bee9396e9b8f
testplot3(PlutoVista)

# ╔═╡ cefb38c1-159e-42db-8088-294573fcece2
md"""
### Changing data

Generally, as above, with Plutovista we need two cells - one with the graph shown, and a second one which triggers the modification.
"""

# ╔═╡ a9f4f98f-ec2f-42d6-88da-4a8a6f727e93
begin
    X = 0:0.05:10
    grid = simplexgrid(X, X)
    f(t) = map((x, y) -> sin(x - t) * atan(y) * cos((y - t)), grid)
end

# ╔═╡ faa59bbd-df1f-4c62-9a77-c4752c6a6df4
vis = GridVisualizer(; resolution = (300, 300), Plotter = PlutoVista, dim = 2);
vis;

# ╔═╡ 6f1707ed-79ab-42dc-8ad8-d66a9e1a65b3
md"""
t= $(@bind t PlutoUI.Slider(-10:0.1:10, default=0, show_value=true))
"""

# ╔═╡ 412c905f-050c-4b78-a66f-0d03978e7edf
scalarplot!(vis, grid, f(t); limits = (-π / 2, π / 2), show = true, levels = 4)

# ╔═╡ e9bc2dae-c303-4063-9ea9-36f95f93371c
md"""
## 2D Grid plot
"""

# ╔═╡ 2b3cb0f4-0656-4981-bec6-48785caf2994
function testgridplot2d(Plotter)
    X = -1:0.2:1
    grid = simplexgrid(X, X)
    return gridplot(grid; Plotter = Plotter, resolution = (300, 300))
end

# ╔═╡ db488643-2e6b-40d6-ba81-126c752953c5
testgridplot2d(PyPlot)

# ╔═╡ 1388c246-be49-4757-a2cc-a685642b6b37
testgridplot2d(PlutoVista)

# ╔═╡ 3f807dfe-7ad3-4740-9e4f-23a39c15603c
testgridplot2d(GLMakie)

# ╔═╡ 5eee8f1d-49ca-4e95-bd14-fe415b0c15e5
md"""
## 3D Scalar plot

Here we use the possibility to update plots to allow moving isosurfaces and plane cuts.
"""

# ╔═╡ 0c99daca-f9a8-4116-867b-e13461c3e754
function grid3d(; n = 15)
    X = collect(0:(1 / n):1)
    return g = simplexgrid(X, X, X)
end

# ╔═╡ 82ccfd24-0053-4399-9bc8-b2e4010bbc92
function func3d(; n = 15)
    g = grid3d(; n = n)
    return g, map((x, y, z) -> sinpi(2 * x) * sinpi(3.5 * y) * sinpi(1.5 * z), g)
end

# ╔═╡ 8b20f720-5470-4da7-bbb6-b746e887046e
g3, f3 = func3d(; n = 19)

# ╔═╡ c0a0ea34-6fc3-4409-934e-086a1a36f94e
p3d = GridVisualizer(; Plotter = PlutoVista, resolution = (300, 300), dim = 3);
p3d;

# ╔═╡ 35be5ef4-0664-4196-8f10-cf71ec7cb371
md"""
f: $(@bind flevel Slider(0:0.01:1,show_value=true,default=0.45))

x: $(@bind xplane Slider(0:0.01:1,show_value=true,default=0.45))
y: $(@bind yplane Slider(0:0.01:1,show_value=true,default=0.45))
z: $(@bind zplane Slider(0:0.01:1,show_value=true,default=0.45))
"""

# ╔═╡ ecd941a0-85b7-4bb7-a903-b19a500198e1
scalarplot!(
    p3d, g3, f3; levels = [flevel], xplanes = [xplane], yplanes = [yplane],
    zplanes = [zplane], colormap = :hot, outlinealpha = 0.05, show = true,
    levelalpha = 0.5
)

# ╔═╡ d924d90d-4102-4ae8-b8de-254a17a5d4df
X4 = -1:0.1:1;
g4 = simplexgrid(X4, X4, X4);

# ╔═╡ 57ed5eea-bc1c-45eb-b4d3-dc63088db21a
scalarplot(
    g4, map((x, y, z) -> 0.01 * exp(-10 * (x^2 + y^2 + z^2)), g4); levels = 3,
    Plotter = PlutoVista
)

# ╔═╡ 597849e9-b9a7-4728-a278-7571d7c1a625
scalarplot(;
    Plotter = PyPlot, g3, f3, resolution = (300, 300), levels = [0.5],
    xplanes = [0.4], yplanes = [0.4], zplanes = [0.4], show = true,
    colormap = :viridis
)

# ╔═╡ 4b9113d2-10bd-4f7a-a2b8-22092656c6b3
md"""
## 3D grid plot
"""

# ╔═╡ 81f0a07d-3d0c-4e7a-9684-1ca4d584b210
gridplot(;
    Plotter = PyPlot, g3, resolution = (300, 300),
    xplanes = [1.0], yplanes = [1.0], zplanes = [0.4], show = true
)

# ╔═╡ f78196ca-d972-4fa6-bdc2-e76eba7ca5a1
p3dgrid = GridVisualizer(; Plotter = PlutoVista, resolution = (300, 300), dim = 3)

# ╔═╡ 7dd92757-c100-4158-baa8-1d9218c39aa7
md"""
x: $(@bind gxplane Slider(0:0.01:1,show_value=true,default=1.0))
y: $(@bind gyplane Slider(0:0.01:1,show_value=true,default=1.0))
z: $(@bind gzplane Slider(0:0.01:1,show_value=true,default=0.45))
"""

# ╔═╡ 840c80b7-5186-45a5-a8df-ec4fb50a5dbb
gridplot!(p3dgrid, g3; xplane = gxplane, yplane = gyplane, zplane = gzplane, show = true)

# ╔═╡ e656dac5-466e-4c07-acfa-0478ad000cb2
function qv2d(; Plotter = nothing, n = 20, stream = false, kwargs...)
    X = 0.0:(10 / (n - 1)):10

    f(x, y) = sin(x) * cos(y) + 0.05 * x * y
    fx(x, y) = -cos(x) * cos(y) - 0.05 * y
    fy(x, y) = sin(x) * sin(y) - 0.05 * x

    grid = simplexgrid(X, X)

    v = map(f, grid)
    ∇v = vcat(map(fx, grid)', map(fy, grid)')

    gvis = GridVisualizer(; Plotter = Plotter, resolution = (600, 600))
    scalarplot!(gvis, grid, v; colormap = :summer, levels = 7)
    if stream
        streamplot!(gvis, grid, ∇v; clear = false, show = true, kwargs...)
    else
        vectorplot!(gvis, grid, ∇v; clear = false, show = true, kwargs...)
    end
    return reveal(gvis)
end

# ╔═╡ 812af347-7606-4c54-b155-88322d20d921
qv2d(; Plotter = PlutoVista, n = 50, spacing = 0.35)

# ╔═╡ 9462f2a4-e7ef-469f-85d4-12f3fa411091
qv2d(; Plotter = PyPlot, n = 100, spacing = 0.35)

# ╔═╡ 8293e0a4-ae15-4929-b20d-95369788f3c9
qv2d(; Plotter = PyPlot, n = 100, spacing = 0.06, stream = true, color = :darkred)

# ╔═╡ 7c66ca9f-7365-489e-ae50-f8c3f9eb145b
qv2d(; Plotter = GLMakie)

# ╔═╡ f03cb172-82c7-4187-be43-01d307977713
qv2d(; Plotter = Plots);

# ╔═╡ Cell order:
# ╠═6df3beed-24a7-4b26-a315-0520f4863190
# ╠═9701cbe0-d048-11eb-151b-67dda7b72b71
# ╠═62c6c85c-e374-4315-b737-b7a0005d2985
# ╠═68e2c958-b417-4ba1-9577-697304fe140a
# ╟─b35d982d-1fa9-413d-b008-892b4f241097
# ╟─00b04f6b-34a6-4e30-9864-d273305281d4
# ╠═a20d74c9-16da-408a-b247-0c17321888f9
# ╠═cc17187f-404c-4c31-8625-fa067eea7273
# ╠═ad1c4dd6-7c6b-4433-a1ac-b2f817ba5d81
# ╠═33482af8-3542-4723-ae43-770a789b69b3
# ╠═d21982eb-aa7e-424c-932e-f0b1460f19b9
# ╠═c4eeb06f-932e-4acc-8e5b-f2a7f9242a42
# ╠═29ca4775-6ba5-474c-bd2c-8f770b09addd
# ╠═c0ab77e8-01ea-436d-85f2-34e253944f11
# ╠═4f73fb99-fab9-4ca4-addb-7d12d16ceba0
# ╠═84192945-d4b6-4949-8f06-d94e04a7a56d
# ╟─7fbaf93f-3cfb-47d0-8252-487e60ba3e54
# ╟─63fe3259-7d79-40ec-98be-e0592e40ee6b
# ╠═4de6b5c9-4d2d-4bcb-bc88-c6f50a23f9b6
# ╟─2061e7fd-c740-4d4b-af5b-7a3a9444aafd
# ╠═f84beb4f-4136-4e5a-ba43-279b703fc75f
# ╟─29fa4467-65ee-4dad-a660-5197864ddbdc
# ╠═c1278fb2-3e75-445f-893a-b8b8a7e931d3
# ╠═661531f7-f740-4dd4-9a59-89ddff06ba5c
# ╠═ed9b80e5-9678-4ba6-bb36-c2e0674ed9ba
# ╠═9ce4f63d-cd96-48d7-a637-07cb84fa88ab
# ╠═77eeefc7-e416-426b-8f87-1bc8439dae6d
# ╠═d503ee1e-1e1f-4235-b286-dc3137a2c96a
# ╠═395d05f4-46a9-4e85-804c-35c79acc96e3
# ╟─ae1fe1ab-4a0e-4c80-bd6f-912201fb4bb4
# ╠═d5258595-60e4-406f-a71e-69111cdad8b9
# ╠═c98a90bf-1a3e-4681-a3b0-663c6844df6b
# ╠═a0c3067b-3aa5-493e-b132-89746483b5ce
# ╠═0998a9a7-d57a-476e-aacd-bee9396e9b8f
# ╟─cefb38c1-159e-42db-8088-294573fcece2
# ╠═a9f4f98f-ec2f-42d6-88da-4a8a6f727e93
# ╠═faa59bbd-df1f-4c62-9a77-c4752c6a6df4
# ╠═412c905f-050c-4b78-a66f-0d03978e7edf
# ╟─6f1707ed-79ab-42dc-8ad8-d66a9e1a65b3
# ╟─e9bc2dae-c303-4063-9ea9-36f95f93371c
# ╠═2b3cb0f4-0656-4981-bec6-48785caf2994
# ╠═db488643-2e6b-40d6-ba81-126c752953c5
# ╠═1388c246-be49-4757-a2cc-a685642b6b37
# ╠═3f807dfe-7ad3-4740-9e4f-23a39c15603c
# ╟─5eee8f1d-49ca-4e95-bd14-fe415b0c15e5
# ╠═0c99daca-f9a8-4116-867b-e13461c3e754
# ╠═82ccfd24-0053-4399-9bc8-b2e4010bbc92
# ╠═8b20f720-5470-4da7-bbb6-b746e887046e
# ╠═c0a0ea34-6fc3-4409-934e-086a1a36f94e
# ╟─35be5ef4-0664-4196-8f10-cf71ec7cb371
# ╠═ecd941a0-85b7-4bb7-a903-b19a500198e1
# ╠═d924d90d-4102-4ae8-b8de-254a17a5d4df
# ╠═57ed5eea-bc1c-45eb-b4d3-dc63088db21a
# ╠═597849e9-b9a7-4728-a278-7571d7c1a625
# ╟─4b9113d2-10bd-4f7a-a2b8-22092656c6b3
# ╠═81f0a07d-3d0c-4e7a-9684-1ca4d584b210
# ╠═f78196ca-d972-4fa6-bdc2-e76eba7ca5a1
# ╟─7dd92757-c100-4158-baa8-1d9218c39aa7
# ╠═840c80b7-5186-45a5-a8df-ec4fb50a5dbb
# ╠═e656dac5-466e-4c07-acfa-0478ad000cb2
# ╠═812af347-7606-4c54-b155-88322d20d921
# ╠═9462f2a4-e7ef-469f-85d4-12f3fa411091
# ╠═8293e0a4-ae15-4929-b20d-95369788f3c9
# ╠═7c66ca9f-7365-489e-ae50-f8c3f9eb145b
# ╠═f03cb172-82c7-4187-be43-01d307977713
