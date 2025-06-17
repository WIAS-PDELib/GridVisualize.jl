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

## Makie
```@docs
GridVisualizeMakieExt
GridVisualizeMakieExt.makescene2d
GridVisualizeMakieExt.makescene3d
GridVisualizeMakieExt.makeaxis3d
GridVisualizeMakieExt.scene_interaction
```
### FlippableLayout

```@docs
GridVisualizeMakieExt.FlippableLayout
GridVisualizeMakieExt.FlippableLayout.FLayout
GridVisualizeMakieExt.FlippableLayout.setmakie!
GridVisualizeMakieExt.FlippableLayout.yieldwait
GridVisualizeMakieExt.FlippableLayout.flayoutscene
```

## PlutoVista
```@docs
GridVisualizePlutoVistaExt
```

## PyPlot
```@docs
GridVisualizePyPlotExt
GridVisualizePyPlotExt.tridata
```
## Plots
```@docs
GridVisualizePlotsExt
GridVisualizePlotsExt.rectdata
```

## ImplEvalSlice
```@docs
GridVisualize.ImplEvalSlice
```

## Experimental
```@docs
GridVisualizeVTKViewExt
GridVisualizeMeshCatExt
```
