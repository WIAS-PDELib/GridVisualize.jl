### A Pluto.jl notebook ###
# v0.20.24

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 9701cbe0-d048-11eb-151b-67dda7b72b71
begin
    using Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "docs"))
    using Revise
    using PlutoVista
    using GridVisualize
    using ExtendableGrids
    using PlutoUI
    using HypertextLiteral
    GridVisualize.default_plotter!(PlutoVista)
end

# ╔═╡ f6205299-d097-4e78-8488-b088475191f6
let
    if isfile("plotting.jl")
        include("plotting.jl")
        plotting_multiscene(; Plotter = PlutoVista, resolution = (650, 300))
    end
end

# ╔═╡ 68e2c958-b417-4ba1-9577-697304fe140a
TableOfContents(; title = "")

# ╔═╡ b35d982d-1fa9-413d-b008-892b4f241097
md"""
# Test notebook for PlutoVista backend of GridVisualize
"""

# ╔═╡ 00b04f6b-34a6-4e30-9864-d273305281d4
md"""
## 1D scalar plot
"""

# ╔═╡ a20d74c9-16da-408a-b247-0c17321888f9
function testplot1()
    grid = simplexgrid(0:0.01:10)
    return scalarplot(
        grid, map(sin, grid); resolution = (600, 200), markershape = :star5,
        markevery = 20, xlabel = "x", ylabel = "z", legend = :rt, label = "sin"
    )
end

# ╔═╡ 33482af8-3542-4723-ae43-770a789b69b3
testplot1()

# ╔═╡ c4eeb06f-932e-4acc-8e5b-f2a7f9242a42
function testplot2(; t = 0)
    p = GridVisualizer(; resolution = (500, 200), legend = :rt, xlabel = "x")
    grid = simplexgrid(0:0.01:10)
    scalarplot!(
        p, grid, map(x -> sin(x - t), grid); color = :red, label = "sin(x-$(t))",
        linestyle = :dash
    )
    scalarplot!(
        p, grid, map(cos, grid); color = :green, clear = false, label = "cos",
        linestyle = :dashdot, linewidth = 3
    )
    return reveal(p)
end

# ╔═╡ 84192945-d4b6-4949-8f06-d94e04a7a56d
testplot2()

# ╔═╡ 7fbaf93f-3cfb-47d0-8252-487e60ba3e54
md"""
### Changing data by re-creating a plot
"""

# ╔═╡ 63fe3259-7d79-40ec-98be-e0592e40ee6b
@bind t2 PlutoUI.Slider(0:0.1:5, show_value = true)

# ╔═╡ 4de6b5c9-4d2d-4bcb-bc88-c6f50a23f9b6
testplot2(; t = t2)

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
p = GridVisualizer(; resolution = (600, 200), dim = 1, legend = :lt)

# ╔═╡ 661531f7-f740-4dd4-9a59-89ddff06ba5c
scalarplot!(p, X2, f2(t4); show = true, clear = true, color = :red, label = "t=$(t4)")

# ╔═╡ dda4599d-05a2-4131-899a-42a653a18b51
md"""
### Non-continuous functions
"""

# ╔═╡ 00b115d3-aa8e-43ef-be6b-3d9d7b42f8af
let
    X = 0:1:10
    g = simplexgrid(X)
    cellmask!(g, [0], [5], 2)
    g1 = subgrid(g, [1])
    g2 = subgrid(g, [2])

    vis = GridVisualizer(; color = :red, size = (600, 200))
    func1 = map((x) -> x, g1)
    func2 = map((x) -> -x, g2)
    func = map(x -> x^2 / 100, g)
    scalarplot!(
        vis,
        [g1, g2],
        g,
        [func1, func2];
        elevation = 0.1,
        clear = true,
        color = :red,
    )
    scalarplot!(
        vis,
        g,
        func;
        elevation = 0.1,
        clear = false,
        color = :green,
    )
    reveal(vis)
end

# ╔═╡ ed9b80e5-9678-4ba6-bb36-c2e0674ed9ba
md"""
## 1D grid plot 
"""

# ╔═╡ 9ce4f63d-cd96-48d7-a637-07cb84fa88ab
function testgridplot()
    grid = simplexgrid(0:0.5:10)
    cellmask!(grid, [0.0], [5], 2)
    bfacemask!(grid, [5], [5], 3)
    return gridplot(grid; resolution = (600, 200), legend = :rt)
end

# ╔═╡ d503ee1e-1e1f-4235-b286-dc3137a2c96a
testgridplot()

# ╔═╡ ae1fe1ab-4a0e-4c80-bd6f-912201fb4bb4
md"""
## 2D Scalar plot
"""

# ╔═╡ d5258595-60e4-406f-a71e-69111cdad8b9
function testplot3()
    X = 0:0.1:10
    grid = simplexgrid(X, X)
    f = map((x, y) -> sin(x) * atan(y), grid)
    return scalarplot(
        grid, f;
        resolution = (300, 300), limits = (-π / 2, π / 2)
    )
end

# ╔═╡ 0998a9a7-d57a-476e-aacd-bee9396e9b8f
testplot3()

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
vis = GridVisualizer(; resolution = (300, 300), dim = 2)

# ╔═╡ 6f1707ed-79ab-42dc-8ad8-d66a9e1a65b3
md"""
t= $(@bind t PlutoUI.Slider(-10:0.1:10, default=0, show_value=true))
"""

# ╔═╡ 412c905f-050c-4b78-a66f-0d03978e7edf
scalarplot!(vis, grid, f(t); limits = (-π / 2, π / 2), show = true, levels = 5)

# ╔═╡ e3c5a486-ec9b-4010-901b-07f7ac997355
md"""
### Non-continuous functions
"""

# ╔═╡ 9b41e77b-e6b9-43b2-89ce-14a9c0eb1242
let
    X = 0:0.7:10
    g = simplexgrid(X, X)
    cellmask!(g, [0, 0], [5, 5], 2)
    g1 = subgrid(g, [1])
    g2 = subgrid(g, [2])

    func1 = map((x, y) -> x^2 + y, g1)
    func2 = map((x, y) -> (x + y^2), g2)
    scalarplot([g1, g2], g, [func1, func2])
end

# ╔═╡ e9bc2dae-c303-4063-9ea9-36f95f93371c
md"""
## 2D Grid plot
"""

# ╔═╡ 2b3cb0f4-0656-4981-bec6-48785caf2994
function testgridplot2d()
    X = -1:0.2:1
    grid = simplexgrid(X, X)
    return gridplot(grid; resolution = (300, 300))
end

# ╔═╡ 1388c246-be49-4757-a2cc-a685642b6b37
testgridplot2d()

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
g3, f3 = func3d(; n = 31)

# ╔═╡ c0a0ea34-6fc3-4409-934e-086a1a36f94e
p3d = GridVisualizer(; resolution = (500, 500), dim = 3)

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
begin
    X4 = -1:0.1:1
    g4 = simplexgrid(X4, X4, X4)
end

# ╔═╡ 57ed5eea-bc1c-45eb-b4d3-dc63088db21a
scalarplot(g4, map((x, y, z) -> 0.01 * exp(-10 * (x^2 + y^2 + z^2)), g4); levels = 10)

# ╔═╡ 943da8f0-d18f-40d5-8158-a3ab5793112f
md"""
### Non-continuous functions
"""

# ╔═╡ ef973737-5cc3-4a3c-8859-a86d9c12c976
let
    X = 0:0.1:1
    g = simplexgrid(X, X, X)
    cellmask!(g, [0, 0, 0], [0.5, 0.5, 0.5], 2)
    g1 = subgrid(g, [1])
    g2 = subgrid(g, [2])
    func1 = map((x, y, z) -> (x + y + z), g1)
    func2 = map((x, y, z) -> (3 - x - y - z), g2)
    scalarplot(
        [g1, g2],
        g,
        [func1, func2];
        levels = 0,
        yplane = 0.25,
        xplane = 0.25,
        zplane = 0.25,
        levelalpha = 1,
        colormap = :hot,
    )
end

# ╔═╡ 4b9113d2-10bd-4f7a-a2b8-22092656c6b3
md"""
## 3D grid plot
"""

# ╔═╡ f78196ca-d972-4fa6-bdc2-e76eba7ca5a1
p3dgrid = GridVisualizer(; resolution = (300, 300), dim = 3)

# ╔═╡ 7dd92757-c100-4158-baa8-1d9218c39aa7
md"""
x: $(@bind gxplane Slider(0:0.01:1,show_value=true,default=1.0))
y: $(@bind gyplane Slider(0:0.01:1,show_value=true,default=1.0))
z: $(@bind gzplane Slider(0:0.01:1,show_value=true,default=0.45))
"""

# ╔═╡ 840c80b7-5186-45a5-a8df-ec4fb50a5dbb
gridplot!(p3dgrid, g3; xplane = gxplane, yplane = gyplane, zplane = gzplane, show = true)

# ╔═╡ 4121e791-8785-472e-a706-7b9cefd36fd6
md"""
## 2D Vector plot
"""

# ╔═╡ e656dac5-466e-4c07-acfa-0478ad000cb2
function qv2d(; n = 20, stream = false, kwargs...)
    X = 0.0:(10 / (n - 1)):10

    f(x, y) = sin(x) * cos(y) + 0.05 * x * y
    fx(x, y) = -cos(x) * cos(y) - 0.05 * y
    fy(x, y) = sin(x) * sin(y) - 0.05 * x

    grid = simplexgrid(X, X)

    v = map(f, grid)
    ∇v = vcat(map(fx, grid)', map(fy, grid)')

    gvis = GridVisualizer(; resolution = (400, 400))
    scalarplot!(gvis, grid, v; colormap = :summer, levels = 7)
    vectorplot!(gvis, grid, ∇v; clear = false, show = true, kwargs...)
    return reveal(gvis)
end

# ╔═╡ 812af347-7606-4c54-b155-88322d20d921
qv2d(; n = 50, rasterpoints = 16)

# ╔═╡ 09998521-68b6-45b4-8c1d-ae73bbd431ad
md"""
## Arranging plots
"""

# ╔═╡ 2d31f310-0d59-4ceb-9daf-61f447de3bb0
let vis = GridVisualizer(; resolution = (600, 600), layout = (2, 2))
    g2 = simplexgrid(X, X)
    g3 = simplexgrid(X, X, X)
    scalarplot!(vis[1, 1], X, sin.(2 * X))
    scalarplot!(vis[1, 2], g2, (x, y) -> sin(x) * cos(y))
    scalarplot!(vis[2, 1], g2, (x, y) -> sin(x) * cos(y); colormap = :hot)
    scalarplot!(
        vis[2, 2], g2, (x, y) -> sin(x) * cos(y); colormap = :hot,
        backend = :plotly
    )
    reveal(vis)
end

# ╔═╡ 6915cc3e-ad9b-4721-9933-884cfc68a25a
md"""
It is also possible to arrange plots using HypertextLiteral:
"""

# ╔═╡ 49db8b25-50ce-4fb4-bea2-de8abfb53c56
X0 = 0:0.25:10

# ╔═╡ 15f4eeb3-c42e-449c-9161-f1df66de6cef
htl"""
<div style= "width: 650px; display: inline-block; white-space:nowrap;">
<div style= "display: inline-block;">
$(scalarplot(X0,X0,X0, (x,y,z)->(sin(x)*sin(y)*sin(z)*sqrt(x*y*z)),resolution=(150,150),colormap=:rainbow))
</div>
<div style= "display: inline-block;">
$(scalarplot(X0,X0,X0, (x,y,z)->(sin(0.1*x*y)*sin(z)),resolution=(150,150),colormap=:hot,backend=:plotly))
</div>
<div style= "display: inline-block;">
$(scalarplot(X0,X0, (x,y)->sin(0.1*x*y),resolution=(150,150),colormap=:summer))
</div>
<div style= "display: inline-block;">
$(scalarplot(X0, (x)->x*sin(2x),color=:red,resolution=(150,150),colormap=:summer))
</div>
</div>
"""

# ╔═╡ 75ffcd09-dfa8-42df-a3cd-a7e68786e73c
md"""
## Misc functionality + tests
"""

# ╔═╡ cf592b99-d596-4511-adbf-001145a59983
md"""
### Plotting of constants
"""

# ╔═╡ bc1c6d12-8d06-4f57-9044-8b5e86fd1c13
scalarplot(ones(10); size = (500, 200))

# ╔═╡ b9c9e4c8-9f4c-481c-bab6-f1baea33c108
md"""
### Aspect ratio handling
"""

# ╔═╡ 6a3b2356-e8a1-45f8-8648-2eca09a7b258
begin
    XX = 0:0.1:1
    YY = 0:0.1:10
end

# ╔═╡ 608a5704-a84c-4c55-8642-ecddb275dc1b
scalarplot(
    XX, YY, (x, y) -> sin(4x) * 10 * y; aspect = 0.1, xlabel = "aaa",
    size = (300, 300)
)

# ╔═╡ 3efbeb11-eaa4-4fc5-bd5f-b3bdb63e7772
scalarplot(
    XX, YY, (x, y) -> sin(4x) * 10 * y; aspect = 0.1, xlabel = "aaa",
    size = (300, 300), backend = :plotly
)

# ╔═╡ ccd274d2-68c0-40e0-8ba7-b8421f5ec9d3
gridplot(simplexgrid(XX, YY); aspect = 0.1)

# ╔═╡ ba5111b8-0dca-42d2-970f-1e88f5392324
html"""<hr>"""

# ╔═╡ Cell order:
# ╠═9701cbe0-d048-11eb-151b-67dda7b72b71
# ╟─68e2c958-b417-4ba1-9577-697304fe140a
# ╟─b35d982d-1fa9-413d-b008-892b4f241097
# ╟─00b04f6b-34a6-4e30-9864-d273305281d4
# ╠═a20d74c9-16da-408a-b247-0c17321888f9
# ╠═33482af8-3542-4723-ae43-770a789b69b3
# ╠═c4eeb06f-932e-4acc-8e5b-f2a7f9242a42
# ╠═84192945-d4b6-4949-8f06-d94e04a7a56d
# ╟─7fbaf93f-3cfb-47d0-8252-487e60ba3e54
# ╟─63fe3259-7d79-40ec-98be-e0592e40ee6b
# ╠═4de6b5c9-4d2d-4bcb-bc88-c6f50a23f9b6
# ╟─2061e7fd-c740-4d4b-af5b-7a3a9444aafd
# ╠═f84beb4f-4136-4e5a-ba43-279b703fc75f
# ╟─29fa4467-65ee-4dad-a660-5197864ddbdc
# ╠═c1278fb2-3e75-445f-893a-b8b8a7e931d3
# ╠═661531f7-f740-4dd4-9a59-89ddff06ba5c
# ╟─dda4599d-05a2-4131-899a-42a653a18b51
# ╠═00b115d3-aa8e-43ef-be6b-3d9d7b42f8af
# ╟─ed9b80e5-9678-4ba6-bb36-c2e0674ed9ba
# ╠═9ce4f63d-cd96-48d7-a637-07cb84fa88ab
# ╠═d503ee1e-1e1f-4235-b286-dc3137a2c96a
# ╟─ae1fe1ab-4a0e-4c80-bd6f-912201fb4bb4
# ╠═d5258595-60e4-406f-a71e-69111cdad8b9
# ╠═0998a9a7-d57a-476e-aacd-bee9396e9b8f
# ╟─cefb38c1-159e-42db-8088-294573fcece2
# ╠═a9f4f98f-ec2f-42d6-88da-4a8a6f727e93
# ╠═faa59bbd-df1f-4c62-9a77-c4752c6a6df4
# ╠═412c905f-050c-4b78-a66f-0d03978e7edf
# ╟─6f1707ed-79ab-42dc-8ad8-d66a9e1a65b3
# ╟─e3c5a486-ec9b-4010-901b-07f7ac997355
# ╠═9b41e77b-e6b9-43b2-89ce-14a9c0eb1242
# ╟─e9bc2dae-c303-4063-9ea9-36f95f93371c
# ╠═2b3cb0f4-0656-4981-bec6-48785caf2994
# ╠═1388c246-be49-4757-a2cc-a685642b6b37
# ╟─5eee8f1d-49ca-4e95-bd14-fe415b0c15e5
# ╠═0c99daca-f9a8-4116-867b-e13461c3e754
# ╠═82ccfd24-0053-4399-9bc8-b2e4010bbc92
# ╠═8b20f720-5470-4da7-bbb6-b746e887046e
# ╠═c0a0ea34-6fc3-4409-934e-086a1a36f94e
# ╟─35be5ef4-0664-4196-8f10-cf71ec7cb371
# ╠═ecd941a0-85b7-4bb7-a903-b19a500198e1
# ╠═d924d90d-4102-4ae8-b8de-254a17a5d4df
# ╠═57ed5eea-bc1c-45eb-b4d3-dc63088db21a
# ╟─943da8f0-d18f-40d5-8158-a3ab5793112f
# ╠═ef973737-5cc3-4a3c-8859-a86d9c12c976
# ╟─4b9113d2-10bd-4f7a-a2b8-22092656c6b3
# ╠═f78196ca-d972-4fa6-bdc2-e76eba7ca5a1
# ╟─7dd92757-c100-4158-baa8-1d9218c39aa7
# ╠═840c80b7-5186-45a5-a8df-ec4fb50a5dbb
# ╟─4121e791-8785-472e-a706-7b9cefd36fd6
# ╠═e656dac5-466e-4c07-acfa-0478ad000cb2
# ╠═812af347-7606-4c54-b155-88322d20d921
# ╟─09998521-68b6-45b4-8c1d-ae73bbd431ad
# ╠═2d31f310-0d59-4ceb-9daf-61f447de3bb0
# ╠═f6205299-d097-4e78-8488-b088475191f6
# ╟─6915cc3e-ad9b-4721-9933-884cfc68a25a
# ╟─49db8b25-50ce-4fb4-bea2-de8abfb53c56
# ╠═15f4eeb3-c42e-449c-9161-f1df66de6cef
# ╟─75ffcd09-dfa8-42df-a3cd-a7e68786e73c
# ╟─cf592b99-d596-4511-adbf-001145a59983
# ╠═bc1c6d12-8d06-4f57-9044-8b5e86fd1c13
# ╟─b9c9e4c8-9f4c-481c-bab6-f1baea33c108
# ╠═6a3b2356-e8a1-45f8-8648-2eca09a7b258
# ╠═608a5704-a84c-4c55-8642-ecddb275dc1b
# ╠═3efbeb11-eaa4-4fc5-bd5f-b3bdb63e7772
# ╠═ccd274d2-68c0-40e0-8ba7-b8421f5ec9d3
# ╟─ba5111b8-0dca-42d2-970f-1e88f5392324
