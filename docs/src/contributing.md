## Hints for contributors
Non-experimental backends are: Makie, PyPlot, PlutoVista, partially Plots. These should be checked before
submitting a pull request.

As it is not easy to install a decent CI for graphical interfaces, visual checking can be performed in the following way.
Create an environment which contains ExtendableGrids, GLMakie, PyPlot, PlutoVista, Plots and PlutoUI, and develops GridVisualize.
A good option is to create a "shared environment" `@GridVisualize` (assuming GridVisualize is worked on in `JULIA_PKG_DEVDIR`):
```
$ julia --project=@GridVisualize
$ julia> ] # Enter Pkg mode
$ (@GridVisualize) pkg> add GLMakie, PyPlot, PlutoVista, Plots, PlutoUI, ExtendableGrids
$ (@GridVisualize) pkg> dev GridVisualize
$ julia> using GridVisualize GLMakie, PyPlot,  Plots
$ julia> include("examples/plotting.jl")
$ julia> plotting_multiscene(Plotter=PyPlot)
$ julia> plotting_multiscene(Plotter=Plots)
$ julia> plotting_multiscene(Plotter=GLMakie)
```
For checking the PlutoVista backend, perform in the same environment 
```
julia> using Pluto
julia> ENV["PLUTO_PROJECT"]=Base.active_project()
julia> Pluto.run(notebook="examples/plutovista.jl")

```

