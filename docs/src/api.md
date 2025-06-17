# Public API

## Visualizer
```@docs
GridVisualizer
GridVisualizer(; Plotter=default_plotter() , kwargs...)
SubVisualizer
reveal
save
```

## Plotters
```@docs
default_plotter
default_plotter!
plottertype
PyPlotType
MakieType
PlotsType
PlutoVistaType
VTKViewType
MeshCatType
```



## Plotting grids
```@docs
gridplot
gridplot!
```


## Plotting scalar data
```@docs
scalarplot
scalarplot!
```

## Plotting vector data
```@docs
vectorplot
vectorplot!
streamplot
streamplot!
```

## Custom plots
```@docs
customplot
customplot!
```
## Plotting TriangulateIO
```@docs
plot_triangulateio
plot_triangulateio!
```

## Keyword Arguments
```@docs
available_kwargs
```

## Supporting methods
```@docs
vectorsample
quiverdata
```

## Creating movies
```@docs
movie
```
