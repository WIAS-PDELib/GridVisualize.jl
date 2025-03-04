# Private API


## Common methods
```@autodocs
Modules = [GridVisualize]
Private = true
Public = false
Pages = ["dispatch.jl","common.jl","slice_plots.jl"]
```

```@docs
ispyplot
isplutovista
isplots
ismakie
isvtkview
```


## PyPlot
```@autodocs
Modules = [GridVisualize]
Private = true
Public = false
Pages = ["pyplot.jl"]
```

## Makie
```@autodocs
Modules = [GridVisualize,FlippableLayout]
Private = true
Public = false
Pages = ["makie.jl", "flippablelayout.jl"]
```

```@docs
FlippableLayout
```

## Plots
```@autodocs
Modules = [GridVisualize]
Private = true
Public = false
Pages = [joinpath("src","plots.jl")] # https://github.com/JuliaDocs/Documenter.jl/issues/2639
```

## VTKView
```@autodocs
Modules = [GridVisualize]
Private = true
Public = false
Pages = ["vtkview.jl"]
```

## Internals
```@docs
GridVisualize.ImplEvalSlice
```
