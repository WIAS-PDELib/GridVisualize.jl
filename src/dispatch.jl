global default_backend = nothing

"""
    default_plotter()

Return default plotter backend. By default, this is set to one of `PlutoVista`,`GLMakie`,`PyPlot` and `Plots`
if one of these packages is imported.
"""
function default_plotter()
    global default_backend
    return default_backend
end

"""
    default_plotter!(Module)

Set plotter module or `nothing` as the default plotter backend.
"""
function default_plotter!(Plotter)
    global default_backend = Plotter
    return default_backend
end

"""
$(SIGNATURES)

Heuristically check if Plotter is VTKView
"""
isvtkview(Plotter) = (typeof(Plotter) == Module) && isdefined(Plotter, :StaticFrame)

"""
$(SIGNATURES)

Heuristically check if Plotter is PyPlot
"""
ispyplot(Plotter) = (typeof(Plotter) == Module) && isdefined(Plotter, :Gcf)

"""
$(SIGNATURES)

Heuristically check if  Plotter is Plots
"""
isplots(Plotter) = (typeof(Plotter) == Module) && isdefined(Plotter, :gr)

"""
$(SIGNATURES)

Heuristically check if Plotter is Makie/WGLMakie
"""
ismakie(Plotter) = (typeof(Plotter) == Module) && isdefined(Plotter, :Makie)

"""
$(SIGNATURES)

Heuristically check if Plotter is MeshCat
"""
ismeshcat(Plotter) = (typeof(Plotter) == Module) && isdefined(Plotter, :Visualizer)

"""
$(SIGNATURES)

Heuristically check if Plotter is PlutoVista
"""
isplutovista(Plotter) = (typeof(Plotter) == Module) && isdefined(Plotter, :PlutoVistaPlot)

"""
$(TYPEDEF)

Abstract type for dispatching on plotter
"""
abstract type AbstractPlotterType end

"""
$(TYPEDEF)

Abstract type for dispatching on plotter
"""
abstract type PyPlotType <: AbstractPlotterType end

"""
$(TYPEDEF)

Abstract type for dispatching on plotter
"""
abstract type MakieType <: AbstractPlotterType end

"""
$(TYPEDEF)

Abstract type for dispatching on plotter
"""
abstract type PlotsType <: AbstractPlotterType end

"""
$(TYPEDEF)

Abstract type for dispatching on plotter. Experimental.
"""
abstract type VTKViewType <: AbstractPlotterType end

"""
$(TYPEDEF)

Abstract type for dispatching on plotter. Experimental.
"""
abstract type MeshCatType <: AbstractPlotterType end

"""
$(TYPEDEF)

Abstract type for dispatching on plotter
"""
abstract type PlutoVistaType <: AbstractPlotterType end

"""
$(SIGNATURES)
    
Heuristically detect type of plotter, returns the corresponding abstract type for plotting.
"""
function plottertype(Plotter::Union{Module, Nothing})
    if ismakie(Plotter)
        return MakieType
    elseif isplots(Plotter)
        return PlotsType
    elseif ispyplot(Plotter)
        return PyPlotType
    elseif isvtkview(Plotter)
        return VTKViewType
    elseif ismeshcat(Plotter)
        return MeshCatType
    elseif isplutovista(Plotter)
        return PlutoVistaType
    end
    return Nothing
end

plottername(::Type{MakieType}) = "Makie"
plottername(::Type{PlotsType}) = "Plots"
plottername(::Type{PyPlotType}) = "PyPlot"
plottername(::Type{PlutoVistaType}) = "PlutoVista"
plottername(::Type{VTKViewType}) = "VTKView"
plottername(::Type{MeshCatType}) = "MeshCat"
plottername(::Type{Nothing}) = "nothing"
plottername(p::Union{Module, Nothing}) = plottertype(p) |> plottername

"""
$(TYPEDEF)

A SubVisualizer is just a dictionary which contains plotting information,
including the type of the plotter and its position in the plot.
"""
const SubVisualizer = Union{Dict{Symbol, Any}, Nothing}

#
# Update subplot context from dict
# and provide some backward compatibilitiy switches.
function _update_context!(ctx::SubVisualizer, kwargs)
    for (k, v) in kwargs
        ctx[Symbol(k)] = v
    end

    haskey(ctx, :resolution) ? ctx[:size] = ctx[:resolution] : nothing
    haskey(ctx, :flevel) ? ctx[:levels] = [ctx[:flevel]] : nothing
    haskey(ctx, :flimits) ? ctx[:limits] = ctx[:flimits] : nothing
    haskey(ctx, :xplane) ? ctx[:xplanes] = [ctx[:xplane]] : nothing
    haskey(ctx, :yplane) ? ctx[:yplanes] = [ctx[:yplane]] : nothing
    haskey(ctx, :zplane) ? ctx[:zplanes] = [ctx[:zplane]] : nothing
    haskey(ctx, :alpha) ? ctx[:outlinealpha] = ctx[:alpha] : nothing

    return ctx
end

"""
$(TYPEDEF)

GridVisualizer struct
"""
struct GridVisualizer
    Plotter::Union{Module, Nothing}
    subplots::Array{SubVisualizer, 2}
    context::SubVisualizer
    function GridVisualizer(
            Plotter::Union{Module, Nothing},
            layout::Tuple,
            default::SubVisualizer
        )
        return new(Plotter, [copy(default) for I in CartesianIndices(layout)], copy(default))
    end
end

"""
$(TYPEDSIGNATURES)

Create a  grid visualizer

Plotter: defaults to `default_plotter()` and can be `PyPlot`, `Plots`, `VTKView`, `Makie` or `PlutoVista´.
This pattern allows  to pass the backend as a module to a plot function without heavy default package dependencies.


Depending on the `layout` keyword argument, a 2D grid of subplots is created.
Further `...plot!` commands then plot into one of these subplots:

```julia
vis=GridVisualizer(Plotter=PyPlot, layout=(2,2)
...plot!(vis[1,2], ...)
```

A `...plot`  command just implicitly creates a plot context:

```julia
gridplot(grid, Plotter=PyPlot) 
```

is equivalent to

```julia
vis=GridVisualizer(Plotter=PyPlot, layout=(1,1))
gridplot!(vis,grid) 
reveal(vis)
```

Please note that the return values of all plot commands are specific to the Plotter.

An interactive mode switch key   for GLMakie (`,`)  and  VTKView (`*`) allows to
toggle between "gallery view" showing all plots at once and "focused view" showing only one plot.


Keyword arguments: see [`available_kwargs`](@ref)

"""
function GridVisualizer(; Plotter::Union{Module, Nothing} = default_plotter(), kwargs...)
    default_ctx = Dict{Symbol, Any}(k => v[1] for (k, v) in default_plot_kwargs())
    _update_context!(default_ctx, kwargs)
    layout = default_ctx[:layout]
    if isnothing(Plotter)
        default_ctx = nothing
    end
    p = GridVisualizer(Plotter, layout, default_ctx)
    if !isnothing(Plotter)
        p.context[:Plotter] = Plotter
        for I in CartesianIndices(layout)
            ctx = p.subplots[I]
            i = Tuple(I)
            ctx[:subplot] = i
            ctx[:iplot] = layout[2] * (i[1] - 1) + i[2]
            ctx[:Plotter] = Plotter
            ctx[:GridVisualizer] = p
        end
        initialize!(p, plottertype(Plotter))
    end
    return p
end

function initialize!(p, ::Type{T}) where {T <: AbstractPlotterType}
    error("Missing implementation of initialize!(p,::Type{$(T)})")
end

function Base.show(io::IO, mime::MIME"text/html", p::GridVisualizer)
    if isplutovista(p.Plotter)
        show(io, mime, p.subplots[1][:figure])
    else
        output = """<code>GridVisualizer(Plotter=$(p.Plotter))</code>"""
        write(io, output)
    end
    return nothing
end

"""
$(SIGNATURES)

Return the layout of a GridVisualizer
"""
Base.size(p::GridVisualizer) = size(p.subplots)

"""
$(SIGNATURES)

Return a SubVisualizer
"""
Base.getindex(p::GridVisualizer, i, j) = p.subplots[i, j]

"""
$(SIGNATURES)

Return the type of a plotter.
"""
plottertype(p::GridVisualizer) = plottertype(p.Plotter)

#
# Default context information with help info.
#
function default_plot_kwargs()
    return OrderedDict{Symbol, Pair{Any, String}}(
        :show => Pair(false, "Show plot immediately"),
        :reveal => Pair(false, "Show plot immediately (same as :show)"),
        :clear => Pair(true, "Clear plot before adding new content"),
        :layout => Pair((1, 1), "Layout of plots in window"),
        :size => Pair((500, 500), "Plot window resolution"),
        :legend => Pair(
            :none,
            "Legend (position): one of [:none, :best, :lt, :ct, :rt, :lc, :rc, :lb, :cb, :rb]"
        ),
        :title => Pair("", "Plot title"),
        :xlabel => Pair("x", "x axis label"),
        :ylabel => Pair("y", "y axis label"),
        :zlabel => Pair("z", "z axis label"),
        :xlimits => Pair((1, -1), "x axis limits"),
        :ylimits => Pair((1, -1), "y axis limits"),
        :zlimits => Pair((1, -1), "z axis limits"),
        :limits => Pair((1, -1), "function limits"),
        :xscale => Pair(:identity, "x axis  scale: one of [:log, :identity]"),
        :yscale => Pair(:identity, "y axis  scale: one of [:log, :identity]"),
        :aspect => Pair(1.0, "XY Aspect ratio modification"),
        :fontsize => Pair(20, "Fontsize of titles. All others are relative to it"),
        :linewidth => Pair(2, "linewidth for isolines or 1D plots"),
        :linestyle => Pair(
            :solid,
            "1D Plot linestyle: one of [:solid, :dash, :dot, :dashdot, :dashdotdot]"
        ),
        :markevery => Pair(5, "1D plot marker stride"),
        :markersize => Pair(5, "1D plot marker size"),
        :markershape => Pair(
            :none,
            "1D plot marker shape: one of [:none, :circle, :star5, :diamond, :hexagon, :cross, :xcross, :utriangle, :dtriangle, :rtriangle, :ltriangle, :pentagon, :+, :x]"
        ),
        :color => Pair(RGB(0.0, 0.0, 0.0), "1D plot line color"),
        :cellwise => Pair(false, "1D plots cellwise; unmaintained and can be slow)"),
        :label => Pair("", "1D plot label"),
        :levels => Pair(7, "array of isolevels or number of isolevels for contour plots"),
        :elevation => Pair(0.0, "2D plot height factor for elevation"),
        :colorlevels => Pair(51, "2D/3D contour plot: number of color levels"),
        :colormap => Pair(
            :viridis,
            "2D/3D contour plot color map (any from [ColorSchemes.jl](https://juliagraphics.github.io/ColorSchemes.jl/stable/basics/#Pre-defined-schemes))"
        ),
        :colorbar => Pair(
            :vertical,
            "2D/3D plot colorbar. One of [:none, :vertical, :horizontal]"
        ),
        :colorbarticks => Pair(
            :default,
            "number of ticks in colorbar (:default sets it equal to levels)"
        ),
        :outlinealpha => Pair(0.05, "3D outline surface alpha value"),
        :levelalpha => Pair(0.25, "3D isolevel alpha"),
        :planealpha => Pair(0.25, "3D plane section alpha"),
        :tetxplane_tol => Pair(0.0, "tolerance for tet-plane intersection in 3D"),
        :rasterpoints => Pair(
            16,
            "Number of quiver points resp. half number of streamplot interpolation points in the maximum extent direction. "
        ),
        :offset => Pair(:default, "Offset of quiver grid"),
        :vscale => Pair(1.0, "Vector field scale for quiver grid"),
        :vconstant => Pair(false, "Set all arrow length constant in vector plot"),
        :vnormalize => Pair(true, "Normalize vector field before scaling"),
        :interior => Pair(true, "3D plot interior of grid"),
        :xplanes => Pair([prevfloat(Inf)], "3D x plane positions or number thereof"),
        :yplanes => Pair([prevfloat(Inf)], "3D y plane positions or number thereof"),
        :zplanes => Pair([prevfloat(Inf)], "3D z plane positions or number thereof"),
        :zoom => Pair(1.0, "Zoom level"),
        :gridscale => Pair(1, "Grid scale factor. Will be applied also to planes, spacing"),
        :cellcoloring => Pair(
            :cellregions,
            "Coloring of cells: one of [:cellregions, :pcolors, :partitions]"
        ),
        :azim => Pair(-60, "3D azimuth angle  (in degrees)"),
        :elev => Pair(30, "3D elevation angle  (in degrees)"),
        :perspectiveness => Pair(
            0.25,
            "3D perspective A number between 0 and 1, where 0 is orthographic, and 1 full perspective"
        ),
        :scene3d => Pair(
            :Axis3,
            "3D plot type of Makie scene. Alternative to `:Axis3` is `:LScene`"
        ),
        :viewmode => Pair(:fit, "Axis3d viewmode for Makie plots. Possible values :fit or :free"),
        :fignumber => Pair(1, "Figure number (PyPlot)"),
        :framepos => Pair(1, "Subplot position in frame (VTKView)"),
        :subplot => Pair((1, 1), "Private: Actual subplot"),
        :backend => Pair(:default, "Backend for PlutoVista plot"),
        :dim => Pair(1, "Data dimension for PlutoVista plot"),
        :regions => Pair(:all, "List of regions to plot"),
        :species => Pair(1, "Number of species to plot or number of species in regions"),
        :spacing => Pair(nothing, "Removed from API"),
        :show_colorbar => Pair(true, "Show color bar next to plots"),
        :slice => Pair(nothing, "Plot a dim-1 slice along a hyperplane expression :(αx ± βy [± γz] ± δ)) or a fixed axis pair, e.g., :x => 3")
    )
end

#
# Print default dict for interpolation into docstrings
#
function _myprint(dict)
    lines_out = IOBuffer()
    for (k, v) in dict
        println(lines_out, "  - `$(k)`: $(v[2]). Default: `$(v[1])`\n")
    end
    return String(take!(lines_out))
end

"""
$(SIGNATURES)

Available kwargs for all methods of this package.

$(_myprint(default_plot_kwargs()))
"""
available_kwargs() = println(_myprint(default_plot_kwargs()))

###################################################################################
function gridplot!(ctx, ::Type{T}, ::Type{Val{D}}, grid) where {T <: AbstractPlotterType, D}
    error("Missing implementation of gridplot!(ctx,::Type{$(T)},  ::Type{Val{$D}}, grid)")
end

"""
$(TYPEDSIGNATURES)

Plot grid into subplot in the visualizer. If `[i,j]` is omitted, `[1,1]` is assumed.

Keyword arguments: see [`available_kwargs`](@ref)
"""
function gridplot!(ctx::SubVisualizer, grid::ExtendableGrid; kwargs...)
    _update_context!(ctx, kwargs)
    return gridplot!(ctx, plottertype(ctx[:Plotter]), Val{dim_space(grid)}, grid)
end

"$(TYPEDSIGNATURES)"
function gridplot!(p::GridVisualizer, grid::ExtendableGrid; kwargs...)
    return gridplot!(p[1, 1], grid; kwargs...)
end

"$(TYPEDSIGNATURES)"
function gridplot!(vis::Union{SubVisualizer, GridVisualizer}, X::AbstractVector; kwargs...)
    return gridplot!(vis, simplexgrid(X); kwargs...)
end

"$(TYPEDSIGNATURES)"
function gridplot!(
        vis::Union{SubVisualizer, GridVisualizer},
        X::AbstractVector,
        Y::AbstractVector; kwargs...
    )

    return gridplot!(vis, simplexgrid(X, Y); kwargs...)
end

"$(TYPEDSIGNATURES)"
function gridplot!(
        vis::Union{SubVisualizer, GridVisualizer},
        X::AbstractVector,
        Y::AbstractVector,
        Z::AbstractVector; kwargs...
    )

    return gridplot!(vis, simplexgrid(X, Y, Z); kwargs...)
end

"$(TYPEDSIGNATURES)"
function gridplot!(
        vis::Union{SubVisualizer, GridVisualizer},
        coord::AbstractMatrix,
        cellnodes::AbstractMatrix,
        kwargs...
    )
    return gridplot!(vis, simplexgrid(coord, cellnodes); kwargs...)
end


"""
$(TYPEDSIGNATURES)

Create grid visualizer and plot grid

Keyword arguments: see [`available_kwargs`](@ref)
"""
function gridplot(grid::ExtendableGrid; Plotter = default_plotter(), kwargs...)
    return gridplot!(GridVisualizer(; Plotter = Plotter, show = true, kwargs...), grid)
end


"$(TYPEDSIGNATURES)"
function gridplot(X::AbstractVector; kwargs...)
    return gridplot(simplexgrid(X); kwargs...)
end

"$(TYPEDSIGNATURES)"
function gridplot(
        X::AbstractVector,
        Y::AbstractVector; kwargs...
    )
    return gridplot(simplexgrid(X, Y); kwargs...)
end

"$(TYPEDSIGNATURES)"
function gridplot(
        X::AbstractVector,
        Y::AbstractVector,
        Z::AbstractVector; kwargs...
    )

    return gridplot(simplexgrid(X, Y, Z); kwargs...)
end

"$(TYPEDSIGNATURES)"
function gridplot(
        coord::AbstractMatrix,
        cellnodes::AbstractMatrix;
        kwargs...
    )
    return gridplot(simplexgrid(coord, cellnodes); kwargs...)
end

###################################################################################
function scalarplot!(ctx, ::Type{T}, ::Type{Val{D}}, grids, parentgrid, funcs) where {T <: AbstractPlotterType, D}
    error("Missing implementation of scalarplot!(ctx,::Type{$(T)},  ::Type{Val{$D}}, grids, parentgrid, funcs)")
end


"""
$(TYPEDSIGNATURES)

Plot node vector on grid as P1 FEM function on the triangulation into subplot in the visualizer. If `[i,j]` is omitted, `[1,1]` is assumed.

If instead of the node vector,  a function is given, it will be evaluated on the grid.

If instead of the grid, coordinate vectors are given, a temporary grid is created.

Keyword arguments: see [`available_kwargs`](@ref)
"""
function scalarplot!(ctx::SubVisualizer, grid::ExtendableGrid, func; kwargs...)
    _update_context!(ctx, kwargs)

    # call a specialized function if the user wants to plot only a dim-1 slice of the data
    if haskey(ctx, :slice) && ctx[:slice] !== nothing
        return slice_plot!(ctx, Val{dim_space(grid)}, grid, func)
    else
        return scalarplot!(ctx, plottertype(ctx[:Plotter]), Val{dim_space(grid)}, [grid], grid, [func])
    end
end

"""
$(TYPEDSIGNATURES)

Plot node vectors on subgrids of parent grid as P1 FEM function on the triangulation into subplot in the visualizer.
If `[i,j]` is omitted, `[1,1]` is assumed.
Keyword arguments: see [`available_kwargs`](@ref)
"""
function scalarplot!(
        ctx::SubVisualizer,
        grids::Vector{ExtendableGrid{Tv, Ti}},
        parentgrid::ExtendableGrid{Tv, Ti},
        funcs::AbstractVector;
        kwargs...,
    ) where {Tv, Ti}
    _update_context!(ctx, kwargs)
    if length(grids) != length(funcs)
        error("number of subgrids: $(length(grids)) and number of functions: $(length(funcs)) not equal")
    end
    return scalarplot!(
        ctx,
        plottertype(ctx[:Plotter]),
        Val{dim_space(parentgrid)},
        grids,
        parentgrid,
        funcs
    )
end

"$(TYPEDSIGNATURES)"
function scalarplot!(p::GridVisualizer, grid::ExtendableGrid, func; kwargs...)
    return scalarplot!(p[1, 1], grid, func; kwargs...)
end

"$(TYPEDSIGNATURES)"
function scalarplot!(ctx::SubVisualizer, grid::ExtendableGrid, func::Function; kwargs...)
    return scalarplot!(ctx, grid, map(func, grid); kwargs...)
end

"$(TYPEDSIGNATURES)"
function scalarplot!(
        p::GridVisualizer,
        grids::Vector{ExtendableGrid{Tv, Ti}},
        parentgrid::ExtendableGrid{Tv, Ti},
        funcs::AbstractVector;
        kwargs...,
    ) where {Tv, Ti}
    return scalarplot!(p[1, 1], grids, parentgrid, funcs; kwargs...)
end

"$(TYPEDSIGNATURES)"
function scalarplot!(ctx::SubVisualizer, func::AbstractVector; kwargs...)
    return scalarplot!(ctx, simplexgrid(1:length(func)), func; kwargs...)
end

"$(TYPEDSIGNATURES)"
function scalarplot!(ctx::SubVisualizer, X::AbstractVector, func; kwargs...)
    return scalarplot!(ctx, simplexgrid(X), func; kwargs...)
end

"$(TYPEDSIGNATURES)"
function scalarplot!(
        ctx::SubVisualizer,
        coord::AbstractMatrix,
        cellnodes::AbstractMatrix,
        func;
        kwargs...,
    )
    return scalarplot!(ctx, simplexgrid(coord, cellnodes), func; kwargs...)
end

"$(TYPEDSIGNATURES)"
function scalarplot!(ctx::GridVisualizer, X::AbstractVector, func; kwargs...)
    return scalarplot!(ctx, simplexgrid(X), func; kwargs...)
end

"$(TYPEDSIGNATURES)"
function scalarplot!(
        ctx::GridVisualizer,
        X::AbstractVector,
        Y::AbstractVector,
        func;
        kwargs...,
    )
    return scalarplot!(ctx, simplexgrid(X, Y), func; kwargs...)
end

"$(TYPEDSIGNATURES)"
function scalarplot!(
        ctx::GridVisualizer,
        X::AbstractVector,
        Y::AbstractVector,
        Z::AbstractVector,
        func;
        kwargs...,
    )
    return scalarplot!(ctx, simplexgrid(X, Y, Z), func; kwargs...)
end

"$(TYPEDSIGNATURES)"
function scalarplot!(
        ctx::GridVisualizer,
        coord::AbstractMatrix,
        cellnodes::AbstractMatrix,
        func;
        kwargs...,
    )
    return scalarplot!(ctx, simplexgrid(coord, cellnodes), func; kwargs...)
end


"""
$(TYPEDSIGNATURES)

Plot node vector on grid as P1 FEM function on the triangulation.

If instead of the node vector,  a function is given, it will be evaluated on the grid.
 
If instead of the grid,  vectors for coordinates are given, a grid is created automatically.

For keyword arguments, see [`available_kwargs`](@ref)
"""
function scalarplot(grid::ExtendableGrid, func; Plotter = default_plotter(), kwargs...)
    return scalarplot!(GridVisualizer(; Plotter = Plotter, kwargs...), grid, func; show = true)
end

function scalarplot(
        grids::Vector{ExtendableGrid{Tv, Ti}},
        parentgrid::ExtendableGrid{Tv, Ti},
        funcs::AbstractVector;
        Plotter = default_plotter(),
        kwargs...,
    ) where {Tv, Ti}
    return scalarplot!(
        GridVisualizer(; Plotter = Plotter, kwargs...),
        grids,
        parentgrid,
        funcs;
        show = true,
    )
end

"$(TYPEDSIGNATURES)"
scalarplot(func::AbstractVector; kwargs...) = scalarplot(1:length(func), func; kwargs...)

"$(TYPEDSIGNATURES)"
scalarplot(X::AbstractVector{T}, func; kwargs...) where {T <: Number} = scalarplot(simplexgrid(X), func; kwargs...)

"$(TYPEDSIGNATURES)"
function scalarplot(X::AbstractVector, Y::AbstractVector, func; kwargs...)
    return scalarplot(simplexgrid(X, Y), func; kwargs...)
end

"$(TYPEDSIGNATURES)"
function scalarplot(
        X::AbstractVector,
        Y::AbstractVector,
        Z::AbstractVector,
        func;
        kwargs...,
    )
    return scalarplot(simplexgrid(X, Y, Z), func; kwargs...)
end

"$(TYPEDSIGNATURES)"
function scalarplot(
        coord::AbstractMatrix,
        cellnodes::AbstractMatrix,
        func;
        kwargs...,
    )
    return scalarplot(simplexgrid(coord, cellnodes), func; kwargs...)
end

###################################################################################
function vectorplot!(ctx, ::Type{T}, ::Type{Val{D}}, grid, func) where {T <: AbstractPlotterType, D}
    error("Missing implementation of vectorplot!(ctx,::Type{$(T)},  ::Type{Val{$D}}, grids, grid,func)")
end


"""
$(TYPEDSIGNATURES)

Plot piecewise linear vector field  as quiver plot.
"""
function vectorplot!(ctx::SubVisualizer, grid::ExtendableGrid, func; kwargs...)
    _update_context!(ctx, kwargs)
    if ctx[:spacing] != nothing
        @warn "`spacing` has been removed from keyword arguments, use `rasterpoints` to control spacing"
    end
    return vectorplot!(ctx, plottertype(ctx[:Plotter]), Val{dim_space(grid)}, grid, func)
end

function vectorplot!(ctx, ptype, ::Type{Val{1}}, grid, func::Matrix)
    return scalarplot!(ctx, ptype, Val{1}, [grid], grid, [func[1, :]])
end
function vectorplot!(ctx, ptype, ::Type{Val{1}}, grid, func::Vector)
    return scalarplot!(ctx, ptype, Val{1}, [grid], grid, [func])
end

"$(TYPEDSIGNATURES)"
function vectorplot!(p::GridVisualizer, grid::ExtendableGrid, func; kwargs...)
    return vectorplot!(p[1, 1], grid, func; kwargs...)
end
"$(TYPEDSIGNATURES)"
function vectorplot!(
        ctx::GridVisualizer,
        X::AbstractVector,
        Y::AbstractVector,
        func;
        kwargs...,
    )
    return vectorplot!(ctx, simplexgrid(X, Y), func; kwargs...)
end
"$(TYPEDSIGNATURES)"
function vectorplot!(
        ctx::GridVisualizer,
        X::AbstractVector,
        Y::AbstractVector,
        Z::AbstractVector,
        func;
        kwargs...,
    )
    return vectorplot!(ctx, simplexgrid(X, Y, Z), func; kwargs...)
end
"$(TYPEDSIGNATURES)"
function vectorplot!(
        ctx::GridVisualizer,
        coord::AbstractMatrix,
        cellnodes::AbstractMatrix,
        func;
        kwargs...,
    )
    return vectorplot!(ctx, simplexgrid(coord, cellnodes), func; kwargs...)
end


"""
$(TYPEDSIGNATURES)

Plot piecewise linear vector field  as quiver plot.
"""
function vectorplot(grid::ExtendableGrid, func; Plotter = default_plotter(), kwargs...)
    return vectorplot!(GridVisualizer(; Plotter = Plotter, kwargs...), grid, func; show = true)
end
"$(TYPEDSIGNATURES)"
function vectorplot(X::AbstractVector, Y::AbstractVector, func; kwargs...)
    return vectorplot(simplexgrid(X, Y), func; kwargs...)
end
"$(TYPEDSIGNATURES)"
function vectorplot(
        X::AbstractVector,
        Y::AbstractVector,
        Z::AbstractVector,
        func;
        kwargs...,
    )
    return vectorplot(simplexgrid(X, Y, Z), func; kwargs...)
end
"$(TYPEDSIGNATURES)"
function vectorplot(
        coord::AbstractMatrix,
        cellnodes::AbstractMatrix,
        func;
        kwargs...,
    )
    return vectorplot(simplexgrid(coord, cellnodes), func; kwargs...)
end

###################################################################################
function streamplot!(ctx, ::Type{T}, ::Type{Val{D}}, grid, func) where {T <: AbstractPlotterType, D}
    error("Missing implementation of streamplot!(ctx,::Type{$(T)},  ::Type{Val{$D}}, grids, grid,func)")
end


"""
$(TYPEDSIGNATURES)

Plot piecewise linear vector field  as stream plot.
(2D pyplot only)
"""
function streamplot!(ctx::SubVisualizer, grid::ExtendableGrid, func; kwargs...)
    _update_context!(ctx, kwargs)
    if ctx[:spacing] != nothing
        @warn "`spacing` has been removed from keyword arguments, use `rasterpoints` to control spacing"
    end
    return streamplot!(ctx, plottertype(ctx[:Plotter]), Val{dim_space(grid)}, grid, func)
end

"$(TYPEDSIGNATURES)"
function streamplot!(p::GridVisualizer, grid::ExtendableGrid, func; kwargs...)
    return streamplot!(p[1, 1], grid, func; kwargs...)
end
"$(TYPEDSIGNATURES)"
function streamplot!(
        ctx::GridVisualizer,
        X::AbstractVector,
        Y::AbstractVector,
        func;
        kwargs...,
    )
    return streamplot!(ctx, simplexgrid(X, Y), func; kwargs...)
end

"$(TYPEDSIGNATURES)"
function streamplot!(
        ctx::GridVisualizer,
        X::AbstractVector,
        Y::AbstractVector,
        Z::AbstractVector,
        func;
        kwargs...,
    )
    return streamplot!(ctx, simplexgrid(X, Y, Z), func; kwargs...)
end

"$(TYPEDSIGNATURES)"
function streamplot!(
        ctx::GridVisualizer,
        coord::AbstractMatrix,
        cellnodes::AbstractMatrix,
        func;
        kwargs...,
    )
    return streamplot!(ctx, simplexgrid(coord, cellnodes), func; kwargs...)
end


"""
$(TYPEDSIGNATURES)

Plot piecewise linear vector field  as stream plot.
(2D pyplot only)
"""
function streamplot(grid::ExtendableGrid, func; Plotter = default_plotter(), kwargs...)
    return streamplot!(GridVisualizer(; Plotter = Plotter, kwargs...), grid, func; show = true)
end
"$(TYPEDSIGNATURES)"
function streamplot(X::AbstractVector, Y::AbstractVector, func; kwargs...)
    return streamplot(simplexgrid(X, Y), func; kwargs...)
end

"$(TYPEDSIGNATURES)"
function streamplot(
        X::AbstractVector,
        Y::AbstractVector,
        Z::AbstractVector,
        func;
        kwargs...,
    )
    return streamplot(simplexgrid(X, Y, Z), func; kwargs...)
end

"$(TYPEDSIGNATURES)"
function streamplot(
        coord::AbstractMatrix,
        cellnodes::AbstractMatrix,
        func;
        kwargs...,
    )
    return streamplot(simplexgrid(coord, cellnodes), func; kwargs...)
end

###################################################################################
"$(TYPEDSIGNATURES)"
function customplot!(ctx::SubVisualizer, func; kwargs...)
    _update_context!(ctx, kwargs)
    if ctx[:spacing] != nothing
        @warn "`spacing` has been removed from keyword arguments, use `rasterpoints` to control spacing"
    end
    return customplot!(ctx, plottertype(ctx[:Plotter]), func)
end

"""
$(TYPEDSIGNATURES)

Variant for do block syntax.
"""
function customplot!(func, ctx::SubVisualizer; kwargs...)
    return customplot!(ctx, func; kwargs...)
end

"$(TYPEDSIGNATURES)"
function customplot!(p::GridVisualizer, func; kwargs...)
    return customplot!(p[1, 1], func; kwargs...)
end

"""
$(TYPEDSIGNATURES)

Variant for do block syntax.
"""
function customplot!(func, p::GridVisualizer; kwargs...)
    return customplot!(p[1, 1], func; kwargs...)
end

"""
$(TYPEDSIGNATURES)

Custom user plot.
"""
function customplot(func; Plotter = default_plotter(), kwargs...)
    return customplot!(GridVisualizer(; Plotter = Plotter, kwargs...), func; show = true)
end

###################################################################################
function reveal(p, ::Type{T}) where {T <: AbstractPlotterType}
    error("Missing implementation of reveal!(p,::Type{$(T)})")
end

"""
$(TYPEDSIGNATURES)

Finish and show plot. Same as setting `:reveal=true` or `:show=true` in last plot statement
for a context.
"""
reveal(visualizer::GridVisualizer) = reveal(visualizer, plottertype(visualizer.Plotter))


###################################################################################
"""
    movie( func, vis::GridVisualizer; file = nothing, format = "gif", kwargs... )

Record a movie with GLMakie backend. MP4 files and gifs can be created.
"""
function movie end

movie(func, visualizer::GridVisualizer, pltype::Any; kwargs...) = nothing

movie(func, visualizer::GridVisualizer; kwargs...) = movie(func, visualizer, plottertype(visualizer.Plotter); kwargs...)

###################################################################################
"""
    plot_triangulateio(;Plotter, 
                       triangulateio; 
                       voronoi=nothing, 
                       circumcircles = false, 
                       kwargs...)

Plot TriangulateIO struct exported by Triangulate.jl. 
"""
function plot_triangulateio end

function plot_triangulateio(
        Plotter,
        triangulateio;
        voronoi = nothing,
        circumcircles = false,
        kwargs...
    )
    vis = GridVisualizer(; Plotter = Plotter, show = true, kwargs...)
    return plot_triangulateio!(vis, triangulateio; voronoi, circumcircles)
end

"$(TYPEDSIGNATURES)"
function plot_triangulateio!(ctx::SubVisualizer, triangulateio; kwargs...)
    _update_context!(ctx, kwargs)
    return plot_triangulateio!(ctx, plottertype(ctx[:Plotter]), triangulateio; kwargs...)
end

"$(TYPEDSIGNATURES)"
function plot_triangulateio!(p::GridVisualizer, triangulateio; kwargs...)
    return plot_triangulateio!(p[1, 1], triangulateio; kwargs...)
end


###################################################################################
function save(fname, p, ::Type{T}) where {T <: AbstractPlotterType}
    error("Missing implementation of save!(fname,p,::Type{$(T)})")
end

function save(fname, scene, PyPlot, ::Type{T}) where {T <: AbstractPlotterType}
    error("Missing implementation of save!(fname,scene, p,::Type{$(T)})")
end

"""
$(TYPEDSIGNATURES)

Save last plotted figure from visualizer to disk.
"""
function save(fname::String, visualizer::GridVisualizer)
    return save(fname, visualizer, plottertype(visualizer.Plotter))
end

"""
$(TYPEDSIGNATURES)

Save scene returned from [`reveal`](@ref), [`scalarplot`](@ref) or [`gridplot`](@ref)  to disk.
"""
function save(fname::String, scene; Plotter = default_plotter())
    return save(fname, scene, Plotter, plottertype(Plotter))
end

#
# Dummy methods to allow Plotter=nothing
#
_update_context!(::Nothing, kwargs) = nothing
Base.copy(::Nothing) = nothing

gridplot!(ctx::Nothing, grid::ExtendableGrid; kwargs...) = nothing
gridplot!(ctx, ::Type{Nothing}, ::Type{Val{1}}, grid) = nothing
gridplot!(ctx, ::Type{Nothing}, ::Type{Val{2}}, grid) = nothing
gridplot!(ctx, ::Type{Nothing}, ::Type{Val{3}}, grid) = nothing

scalarplot!(ctx::Nothing, grid::ExtendableGrid, func; kwargs...) = nothing
scalarplot!(ctx::Nothing, grid::ExtendableGrid, func::Function; kwargs...) = nothing
scalarplot!(ctx, ::Type{Nothing}, ::Type{Val{1}}, grid, func) = nothing
scalarplot!(ctx, ::Type{Nothing}, ::Type{Val{2}}, grid, func) = nothing
scalarplot!(ctx, ::Type{Nothing}, ::Type{Val{3}}, grid, func) = nothing

vectorplot!(ctx::Nothing, grid::ExtendableGrid, func; kwargs...) = nothing
vectorplot!(ctx::Nothing, grid::ExtendableGrid, func::Function; kwargs...) = nothing
vectorplot!(ctx, ::Type{Nothing}, ::Type{Val{2}}, grid, func) = nothing
vectorplot!(ctx, ::Type{Nothing}, ::Type{Val{3}}, grid, func) = nothing

streamplot!(ctx::Nothing, grid::ExtendableGrid, func; kwargs...) = nothing
streamplot!(ctx::Nothing, grid::ExtendableGrid, func::Function; kwargs...) = nothing
streamplot!(ctx, ::Type{Nothing}, ::Type{Val{2}}, grid, func) = nothing
streamplot!(ctx, ::Type{Nothing}, ::Type{Val{3}}, grid, func) = nothing

customplot!(ctx::Nothing, func::Function; kwargs...) = nothing

save(fname, scene, Plotter::Nothing, ::Type{Nothing}) = nothing
displayable(ctx, Any) = nothing
reveal(p, ::Type{Nothing}) = nothing
