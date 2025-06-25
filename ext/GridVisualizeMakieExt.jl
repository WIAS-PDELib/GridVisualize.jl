"""
    module GridVisualizeMakieExt

Makie extension module for GridVisualize.jl. This is almost fully supported and works with GLMakie,
CairoMalie and WGLMakie.
"""
module GridVisualizeMakieExt

using Colors
using ColorSchemes
using DocStringExtensions
using Printf
using GeometryBasics
using Interpolations: linear_interpolation
using IntervalSets

import GridVisualize: initialize!, save, reveal, gridplot!, scalarplot!, vectorplot!, streamplot!, customplot!, movie, plot_triangulateio!
using GridVisualize: MakieType, GridVisualizer, SubVisualizer
using GridVisualize: isolevels, cellcolors, num_cellcolors, vectorsample, quiverdata, regionmesh, bfacesegments

using ExtendableGrids
using GridVisualizeTools
using Observables

include("flippablelayout.jl")

function initialize!(p::GridVisualizer, ::Type{MakieType})
    XMakie = p.context[:Plotter]

    # Prepare flippable layout
    FlippableLayout.setmakie!(XMakie.Makie)

    layout = p.context[:layout]

    parent, flayout = FlippableLayout.flayoutscene(; size = p.context[:size])

    p.context[:figure] = parent
    p.context[:flayout] = flayout

    # copy arguments to sublayout
    for I in CartesianIndices(layout)
        ctx = p.subplots[I]
        ctx[:figure] = parent
        ctx[:flayout] = flayout
    end

    # Don't call display on pluto
    if !isdefined(Main, :PlutoRunner)
        XMakie.display(parent)
    end

    return parent
end

# Adding a scene to the layout just adds to the
# flippable layout.
add_scene!(ctx, ax) = ctx[:flayout][ctx[:subplot]...] = ax

# Revealing the  visualizer just returns the figure
function reveal(p::GridVisualizer, ::Type{MakieType})
    XMakie = p.context[:Plotter]
    layout = p.context[:layout]
    # For 1D plots the legend should be rendered only once,
    # when all lines+labels are defined
    for I in CartesianIndices(layout)
        ctx = p.subplots[I]
        if ctx[:legend] != :none && haskey(ctx, :scalarplot1d)
            if !haskey(ctx, :axislegend)
                pos = ctx[:legend] == :best ? :rt : ctx[:legend]
                ctx[:axislegend] = XMakie.axislegend(
                    ctx[:scene];
                    position = pos,
                    labelsize = 0.5 * ctx[:fontsize],
                    backgroundcolor = RGBA(1.0, 1.0, 1.0, 0.85),
                )
            end
        end
    end

    if haskey(p.context, :videostream)
        return XMakie.recordframe!(p.context[:videostream])
    else
        return p.context[:figure]
    end
end

function reveal(ctx::SubVisualizer, TP::Type{MakieType})
    FlippableLayout.yieldwait(ctx[:flayout])
    if ctx[:show] || ctx[:reveal]
        return reveal(ctx[:GridVisualizer], TP)
    end
    return nothing
end

function save(fname, p::GridVisualizer, ::Type{MakieType})
    return p.context[:Plotter].save(fname, p.context[:figure])
end
function save(fname, scene, XMakie, ::Type{MakieType})
    return isnothing(scene) ? nothing : XMakie.save(fname, scene)
end


function movie(
        func,
        vis::GridVisualizer,
        ::Type{MakieType};
        file = nothing,
        format = "gif",
        kwargs...,
    )
    Plotter = vis.context[:Plotter]
    if !isnothing(file)
        format = lstrip(splitext(file)[2], '.')
    end

    if isdefined(Main, :PlutoRunner) || !isnothing(file)
        vis.context[:videostream] = Plotter.VideoStream(vis.context[:figure]; format = format, kwargs...)
    end

    func(vis)

    if !isnothing(file)
        return Plotter.save(file, vis.context[:videostream])
    elseif isdefined(Main, :PlutoRunner)
        return vis.context[:videostream]
    else
        return nothing
    end
end

"""

     scene_interaction(update_scene,view,switchkeys::Vector{Symbol}=[:nothing])   

Control multiple scene elements via keyboard up/down keys. 
Each switchkey is assumed to correspond to one of these elements.
Pressing a switch key transfers control to its associated element.

Control of values of the current associated element is performed
by triggering change values via up/down (± 1)  resp. page_up/page_down (±10) keys

The update_scene callback gets passed the change value and the symbol.
"""
function scene_interaction(
        update_scene,
        scene,
        XMakie,
        switchkeys::Vector{Symbol} = [:nothing]
    )

    # Check if pixel position pos sits within the scene
    function _inscene(scene, pos)
        area = scene.viewport[]
        return pos[1] > area.origin[1] &&
            pos[1] < area.origin[1] + area.widths[1] &&
            pos[2] > area.origin[2] &&
            pos[2] < area.origin[2] + area.widths[2]
    end

    # Initial active switch key is the first in the vector passed
    activeswitch = Observable(switchkeys[1])

    # Handle mouse position within scene
    mouseposition = Observable((0.0, 0.0))

    XMakie.on(scene.events.mouseposition) do m
        mouseposition[] = m
        false
    end

    # Set keyboard event callback
    return XMakie.on(scene.events.keyboardbutton) do buttons
        if _inscene(scene, mouseposition[])
            # On pressing a switch key, pass control
            for i in 1:length(switchkeys)
                if switchkeys[i] != :nothing &&
                        XMakie.ispressed(scene, getproperty(XMakie.Keyboard, switchkeys[i]))
                    activeswitch[] = switchkeys[i]
                    update_scene(0, switchkeys[i])
                    return true
                end
            end

            # Handle change values via up/down control
            if XMakie.ispressed(scene, XMakie.Keyboard.up)
                update_scene(1, activeswitch[])
                return true
            elseif XMakie.ispressed(scene, XMakie.Keyboard.down)
                update_scene(-1, activeswitch[])
                return true
            elseif XMakie.ispressed(scene, XMakie.Keyboard.page_up)
                update_scene(10, activeswitch[])
                return true
            elseif XMakie.ispressed(scene, XMakie.Keyboard.page_down)
                update_scene(-10, activeswitch[])
                return true
            end
        end
        return false
    end
end

# Standard kwargs for Makie scenes
scenekwargs(ctx) = Dict(
    #:xticklabelsize => 0.5*ctx[:fontsize],
    #:yticklabelsize => 0.5*ctx[:fontsize],
    #:zticklabelsize => 0.5*ctx[:fontsize],
    #:xlabelsize => 0.5*ctx[:fontsize],
    #:ylabelsize => 0.5*ctx[:fontsize],
    #:zlabelsize => 0.5*ctx[:fontsize],
    #:xlabeloffset => 20,
    #:ylabeloffset => 20,
    #:zlabeloffset => 20,
    :titlesize => ctx[:fontsize]
)

#scenekwargs(ctx)=()

############################################################################################################
#1D grid

# Point list for node markers
function basemesh1d(grid, gridscale)
    coord = vec(grid[Coordinates])
    ncoord = length(coord)
    points = Vector{Point2f}(undef, 0)
    (xmin, xmax) = extrema(coord)
    h = gridscale * (xmax - xmin) / 40.0
    for i in 1:ncoord
        push!(points, Point2f(coord[i] * gridscale, h))
        push!(points, Point2f(coord[i] * gridscale, -h))
    end
    return points
end

# Point list for intervals
function regionmesh1d(grid, gridscale, iregion; cellcoloring = :cellregions)
    coord = vec(grid[Coordinates])
    points = Vector{Point2f}(undef, 0)
    cn = grid[CellNodes]
    cr = cellcolors(grid, cellcoloring)
    ncells = length(cr)
    for i in 1:ncells
        if cr[i] == iregion
            push!(points, Point2f(coord[cn[1, i]] * gridscale, 0))
            push!(points, Point2f(coord[cn[2, i]] * gridscale, 0))
        end
    end
    return points
end

# Point list for boundary nodes
function bregionmesh1d(grid, gridscale, ibreg)
    nbfaces = num_bfaces(grid)
    bfacenodes = grid[BFaceNodes]
    bfaceregions = grid[BFaceRegions]
    coord = vec(grid[Coordinates])
    points = Vector{Point2f}(undef, 0)
    (xmin, xmax) = extrema(coord)
    h = gridscale * (xmax - xmin) / 20.0
    for ibface in 1:nbfaces
        if bfaceregions[ibface] == ibreg
            push!(points, Point2f(coord[bfacenodes[1, ibface]] * gridscale, h))
            push!(points, Point2f(coord[bfacenodes[1, ibface]] * gridscale, -h))
        end
    end
    return points
end

# Point list for scene size
function scenecorners1d(grid, gridscale)
    coord = vec(grid[Coordinates])
    (xmin, xmax) = extrema(coord)
    h = gridscale * (xmax - xmin) / 40.0
    return [Point2f(xmin * gridscale, -5 * h), Point2f(xmax * gridscale, 5 * h)]
end

function gridplot!(ctx, TP::Type{MakieType}, ::Type{Val{1}}, grid)
    XMakie = ctx[:Plotter]
    nregions = num_cellregions(grid)
    nbregions = num_bfaceregions(grid)
    gridscale = ctx[:gridscale]
    if !haskey(ctx, :scene)
        ctx[:scene] = XMakie.Axis(
            ctx[:figure];
            yticklabelsvisible = false,
            yticksvisible = false,
            title = ctx[:title],
            scenekwargs(ctx)...,
        )

        ctx[:grid] = Observable(grid)
        cmap = region_cmap(nregions)
        bcmap = bregion_cmap(nbregions)

        # Set scene size with invisible markers
        XMakie.scatter!(
            ctx[:scene],
            map(g -> scenecorners1d(grid, gridscale), ctx[:grid]);
            color = :white,
            markersize = 0.0,
            strokewidth = 0,
        )

        # Draw node markers
        XMakie.linesegments!(
            ctx[:scene],
            map(g -> basemesh1d(g, gridscale), ctx[:grid]);
            color = :black,
        )

        # Colored cell regions
        for i in 1:nregions
            XMakie.linesegments!(
                ctx[:scene],
                map(g -> regionmesh1d(g, gridscale, i; cellcoloring = ctx[:cellcoloring]), ctx[:grid]);
                color = cmap[i],
                linewidth = 4,
                label = "c $(i)",
            )
        end

        # Colored boundary grid
        for i in 1:nbregions
            XMakie.linesegments!(
                ctx[:scene],
                map(g -> bregionmesh1d(g, gridscale, i), ctx[:grid]);
                color = bcmap[i],
                linewidth = 4,
                label = "b$(i)",
            )
        end

        # Legende
        if ctx[:legend] != :none
            pos = ctx[:legend] == :best ? :rt : ctx[:legend]
            XMakie.axislegend(
                ctx[:scene];
                position = pos,
                labelsize = 0.5 * ctx[:fontsize],
                nbanks = 5,
            )
        end
        XMakie.reset_limits!(ctx[:scene])
        add_scene!(ctx, ctx[:scene])
    else
        ctx[:grid][] = grid
    end
    return reveal(ctx, TP)
end

########################################################################
# 1D function

function scalarplot!(ctx, TP::Type{MakieType}, ::Type{Val{1}}, grids, parentgrid, funcs)
    XMakie = ctx[:Plotter]

    nfuncs = length(funcs)
    if ctx[:title] == ""
        ctx[:title] = " "
    end
    ctx[:scalarplot1d] = true
    gridscale = ctx[:gridscale]
    # ... keep this for the case we are unsorted
    function polysegs(grid, func)
        points = Vector{Point2f}(undef, 0)
        cellnodes = grid[CellNodes]
        coord = grid[Coordinates]
        for icell in 1:num_cells(grid)
            i1 = cellnodes[1, icell]
            i2 = cellnodes[2, icell]
            x1 = coord[1, i1] * gridscale
            x2 = coord[1, i2] * gridscale
            push!(points, Point2f(x1, func[i1]))
            push!(points, Point2f(x2, func[i2]))
        end
        return points
    end

    function polyline(grid, func)
        coord = grid[Coordinates]
        return points = [Point2f(coord[1, i] * gridscale, func[i]) for i in 1:length(func)]
    end

    coord = parentgrid[Coordinates]
    xlimits = ctx[:xlimits]
    ylimits = ctx[:limits]

    xmin = coord[1, 1] * gridscale
    xmax = coord[1, end] * gridscale
    xauto = true
    yauto = true
    if xlimits[1] < xlimits[2]
        xmin = xlimits[1]
        xmax = xlimits[2]
        xauto = false
    end

    if ylimits[1] < ylimits[2]
        ymin = ylimits[1]
        ymax = ylimits[2]
        yauto = false
    else
        ext = extrema.(funcs)
        (ymin, ymax) = (minimum(first.(ext)), maximum(last.(ext)))
    end

    function update_lines(ctx, newrange)
        return if ctx[:markershape] == :none
            #line without marker
            for l in newrange
                XMakie.lines!(
                    ctx[:scene],
                    map(a -> a, ctx[:lines][l]);
                    linestyle = ctx[:linestyle],
                    linewidth = ctx[:linewidth],
                    color = rgbcolor(ctx[:color])
                )
            end
            if ctx[:label] != ""
                XMakie.scatterlines!(
                    ctx[:scene],
                    map(a -> a[1:1], ctx[:lines][newrange[begin]]);
                    linestyle = ctx[:linestyle],
                    linewidth = ctx[:linewidth],
                    markersize = 0.1,
                    color = rgbcolor(ctx[:color]),
                    label = ctx[:label],
                )
            end
        else
            # line with markers separated by markevery

            # draw plain line without the label
            for l in newrange
                XMakie.lines!(
                    ctx[:scene],
                    map(a -> a, ctx[:lines][l]);
                    linestyle = ctx[:linestyle],
                    color = rgbcolor(ctx[:color]),
                    linewidth = ctx[:linewidth],
                )
                # draw markers without label
                XMakie.scatter!(
                    ctx[:scene],
                    map(a -> a[1:ctx[:markevery]:end], ctx[:lines][l]);
                    color = rgbcolor(ctx[:color]),
                    marker = ctx[:markershape],
                    markersize = ctx[:markersize],
                )
            end

            # Draw  dummy line with marker on top of the first
            # marker position already drawn in order to
            # get the proper legend entry
            if ctx[:label] != ""
                XMakie.scatterlines!(
                    ctx[:scene],
                    map(a -> a[1:1], ctx[:lines][newrange[begin]]);
                    linestyle = ctx[:linestyle],
                    linewidth = ctx[:linewidth],
                    marker = ctx[:markershape],
                    markersize = ctx[:markersize],
                    markercolor = rgbcolor(ctx[:color]),
                    color = rgbcolor(ctx[:color]),
                    label = ctx[:label],
                )
            end
        end
    end

    if !haskey(ctx, :scene)
        ctx[:xtitle] = Observable(ctx[:title])

        # Axis
        ctx[:scene] = XMakie.Axis(
            ctx[:figure];
            title = map(a -> a, ctx[:xtitle]),
            xscale = ctx[:xscale] == :log ? log10 : identity,
            yscale = ctx[:yscale] == :log ? log10 : identity,
            xlabel = ctx[:xlabel],
            ylabel = ctx[:ylabel],
            scenekwargs(ctx)...,
        )

        if !xauto
            XMakie.xlims!(ctx[:scene], xmin, xmax)
        end
        if !yauto
            XMakie.ylims!(ctx[:scene], ymin, ymax)
        end
        # Plot size
        XMakie.scatter!(
            ctx[:scene],
            [Point2f(xmin, ymin), Point2f(xmax, ymax)];
            color = parse(RGB{Float64}, :white),
            markersize = 0.0,
            strokewidth = 0,
        )

        # ctx[:lines]  is an array of lines to draw
        # Here, we start just with the first one.

        ctx[:lines] = [Observable(polyline(grids[i], funcs[i])) for i in 1:nfuncs]

        update_lines(ctx, 1:nfuncs)

        XMakie.reset_limits!(ctx[:scene])

        ctx[:nlines] = nfuncs

        XMakie.reset_limits!(ctx[:scene])
        add_scene!(ctx, ctx[:scene])

    else
        if ctx[:clear]
            ctx[:nlines] = nfuncs
        else
            ctx[:nlines] += nfuncs
        end

        # Either update existing line, or
        # create new one. This works with repeating sequences of
        # updating lines.
        if ctx[:nlines] <= length(ctx[:lines])
            for i in 1:nfuncs
                ctx[:lines][ctx[:nlines] - nfuncs + i][] = polyline(grids[i], funcs[i])
            end
        else
            r0 = length(ctx[:lines])
            for i in 1:nfuncs
                push!(ctx[:lines], Observable(polyline(grids[i], funcs[i])))
            end
            r1 = length(ctx[:lines])
            update_lines(ctx, (r0 + 1):r1)
        end

        XMakie.reset_limits!(ctx[:scene])

        ctx[:xtitle][] = ctx[:title]
    end

    return reveal(ctx, TP)
end

#######################################################################################
# 2D grid

function makescene_grid(ctx)
    XMakie = ctx[:Plotter]
    GL = XMakie.GridLayout(ctx[:figure])
    GL[1, 1] = ctx[:scene]
    ncol = length(ctx[:cmap])
    nbcol = length(ctx[:bcmap])
    # fontsize=0.5*ctx[:fontsize],ticklabelsize=0.5*ctx[:fontsize]
    if ctx[:show_colorbar]
        if ctx[:colorbar] == :vertical
            GL[1, 2] = XMakie.Colorbar(
                ctx[:figure];
                colormap = XMakie.cgrad(ctx[:cmap]; categorical = true),
                limits = (0.5, ncol + 0.5),
                ticks = 1:ncol,
                width = 15,
                label = "cell regions",
            )
            GL[1, 3] = XMakie.Colorbar(
                ctx[:figure];
                colormap = XMakie.cgrad(ctx[:bcmap]; categorical = true),
                limits = (0.5, nbcol + 0.5),
                ticks = 1:nbcol,
                width = 15,
                label = "boundary regions",
            )
        elseif ctx[:colorbar] == :horizontal
            GL[2, 1] = XMakie.Colorbar(
                ctx[:figure];
                colormap = XMakie.cgrad(ctx[:cmap]; categorical = true),
                limits = (1, ncol),
                height = 15,
                vertical = false,
                label = "cell regions",
            )
            GL[3, 1] = XMakie.Colorbar(
                ctx[:figure];
                colormap = XMakie.cgrad(ctx[:bcmap]; categorical = true),
                limits = (1, nbcol),
                height = 15,
                vertical = false,
                label = "boundary regions",
            )
        end
    end
    return GL
end

# Put all data which could be updated in to one plot.
function set_plot_data!(ctx, key, data)
    return haskey(ctx, key) ? ctx[key][] = data : ctx[key] = Observable(data)
end

function gridplot!(ctx, TP::Type{MakieType}, ::Type{Val{2}}, grid)
    XMakie = ctx[:Plotter]

    nregions = num_cellcolors(grid, ctx[:cellcoloring])

    nbregions = num_bfaceregions(grid)

    set_plot_data!(ctx, :grid, grid)

    if !haskey(ctx, :gridplot)
        if !haskey(ctx, :scene)
            aspect = nothing
            autolimitaspect = nothing
            if ctx[:aspect] ≈ 1.0
                aspect = XMakie.DataAspect()
            else
                autolimitaspect = ctx[:aspect]
            end
            ctx[:scene] = XMakie.Axis(
                ctx[:figure];
                title = ctx[:title],
                aspect = aspect,
                autolimitaspect = autolimitaspect,
                scenekwargs(ctx)...,
            )
            xlimits = ctx[:xlimits]
            ylimits = ctx[:ylimits]
            if xlimits[1] < xlimits[2]
                XMakie.xlims!(ctx[:scene], xlimits...)
            end
            if ylimits[1] < ylimits[2]
                XMakie.ylims!(ctx[:scene], ylimits...)
            end
        end
        # Draw cells with region mark
        cmap = region_cmap(nregions)
        ctx[:cmap] = cmap
        for i in 1:nregions
            XMakie.poly!(
                ctx[:scene],
                map(g -> regionmesh(g, ctx[:gridscale], i; cellcoloring = ctx[:cellcoloring]), ctx[:grid]);
                color = cmap[i],
                strokecolor = :black,
                strokewidth = ctx[:linewidth],
            )
        end

        # Draw boundary lines
        bcmap = bregion_cmap(nbregions)
        ctx[:bcmap] = bcmap
        for i in 1:nbregions
            lp = XMakie.linesegments!(
                ctx[:scene],
                map(g -> bfacesegments(g, ctx[:gridscale], i), ctx[:grid]);
                label = "$(i)",
                color = bcmap[i],
                linewidth = 4,
            )
            XMakie.translate!(lp, 0, 0, 0.1)
        end
        XMakie.reset_limits!(ctx[:scene])

        # Describe legend
        if ctx[:legend] != :none
            pos = ctx[:legend] == :best ? :rt : ctx[:legend]
            XMakie.axislegend(
                ctx[:scene];
                position = pos,
                labelsize = 0.5 * ctx[:fontsize],
                backgroundcolor = :transparent,
            )
        end
        add_scene!(ctx, makescene_grid(ctx))
    end
    return reveal(ctx, TP)
end

"""
   makescene2d(ctx)

Complete scene with title and status line showing interaction state.
This uses a gridlayout and its  protrusion capabilities.
"""
function makescene2d(ctx, key)
    XMakie = ctx[:Plotter]
    GL = XMakie.GridLayout(ctx[:figure])
    GL[1, 1] = ctx[:scene]

    # , fontsize=0.5*ctx[:fontsize],ticklabelsize=0.5*ctx[:fontsize]
    if ctx[:show_colorbar]
        if ctx[:colorbar] == :vertical
            GL[1, 2] = XMakie.Colorbar(
                ctx[:figure],
                ctx[key];
                width = 10,
                ticks = unique(ctx[:cbarticks]),
                tickformat = "{:.2e}",
            )
        elseif ctx[:colorbar] == :horizontal
            GL[2, 1] = XMakie.Colorbar(
                ctx[:figure],
                ctx[key];
                height = 10,
                ticks = unique(ctx[:cbarticks]),
                vertical = false,
                tickformat = "{:.2e}",
            )
        end
    end
    return GL
end

# 2D function
function scalarplot!(ctx, TP::Type{MakieType}, ::Type{Val{2}}, grids, parentgrid, funcs)
    XMakie = ctx[:Plotter]
    gridscale = ctx[:gridscale]
    regions = ctx[:regions]

    # set limits according to extrema of funcs on the specified grid regions
    if ctx[:limits] == :fit && regions != :all
        f_min, f_max = (Inf64, -Inf64)
        for (func, grid) in zip(funcs, grids)
            for icell in 1:num_cells(grid)
                if grid[CellRegions][icell] in regions
                    @views f_min = min(f_min, func[grid[CellNodes][:, icell]]...)
                    @views f_max = max(f_max, func[grid[CellNodes][:, icell]]...)
                end
            end
        end
        ctx[:limits] = (f_min, f_max)
    end

    # Create GeometryBasics.mesh from grid data.
    function make_mesh(grids, funcs, elevation)
        ngrids = length(grids)
        coords = [grid[Coordinates] for grid in grids]
        npoints = [num_nodes(grid) for grid in grids]
        cellnodes = [grid[CellNodes] for grid in grids]
        ncells = [num_cells(grid) for grid in grids]
        cellregions = [grid[CellRegions] for grid in grids]
        offsets = zeros(Int, ngrids)
        for i in 2:ngrids
            offsets[i] = offsets[i - 1] + npoints[i - 1]
        end

        if elevation ≈ 0.0
            points = Vector{Point3f}(undef, sum(npoints))
            k = 1
            for j in 1:ngrids
                for i in 1:npoints[j]
                    points[k] = Point3f(coords[j][1, i] * gridscale, coords[j][2, i] * gridscale, -0.1)
                    k = k + 1
                end
            end
        else
            points = Vector{Point3f}(undef, sum(npoints))
            k = 1
            for j in 1:ngrids
                for i in 1:npoints[j]
                    points[k] = Point3f(
                        coords[j][1, i] * gridscale, coords[j][2, i] * gridscale,
                        funcs[j][i] * elevation * gridscale
                    )
                    k = k + 1
                end
            end
        end
        faces = Vector{TriangleFace{Int64}}(undef, 0)
        for j in 1:ngrids
            for i in 1:ncells[j]
                if regions == :all || (cellregions[j][i] in regions)
                    push!(
                        faces,
                        TriangleFace(
                            cellnodes[j][1, i] + offsets[j],
                            cellnodes[j][2, i] + offsets[j],
                            cellnodes[j][3, i] + offsets[j]
                        )
                    )
                end
            end
        end
        return Mesh(points, faces)
    end

    levels, crange, ctx[:cbarticks] = isolevels(ctx, funcs)

    eps = 1.0e-1
    if crange[1] == crange[2]
        crange = (crange[1] - eps, crange[1] + eps)
    end

    set_plot_data!(
        ctx,
        :contourdata,
        (
            g = grids,
            f = funcs,
            e = ctx[:elevation],
            t = ctx[:title],
            l = levels,
            c = crange,
        )
    )

    if !haskey(ctx, :contourplot)
        if !haskey(ctx, :scene)
            aspect = nothing
            autolimitaspect = nothing
            if ctx[:aspect] ≈ 1.0
                aspect = XMakie.DataAspect()
            else
                autolimitaspect = ctx[:aspect]
            end
            if ctx[:elevation] ≈ 0
                ctx[:scene] = XMakie.Axis(
                    ctx[:figure];
                    title = map(data -> data.t, ctx[:contourdata]),
                    aspect = aspect,
                    autolimitaspect = autolimitaspect,
                    xlabel = ctx[:xlabel],
                    ylabel = ctx[:ylabel],
                    scenekwargs(ctx)...,
                )
            else
                ctx[:scene] = XMakie.Axis3(
                    ctx[:figure];
                    title = map(data -> data.t, ctx[:contourdata]),
                    aspect = aspect,
                    autolimitaspect = autolimitaspect,
                    scenekwargs(ctx)...,
                )
            end

            xlimits = ctx[:xlimits]
            ylimits = ctx[:ylimits]
            if xlimits[1] < xlimits[2]
                XMakie.xlims!(ctx[:scene], xlimits...)
            end
            if ylimits[1] < ylimits[2]
                XMakie.ylims!(ctx[:scene], ylimits...)
            end
        end

        # Draw the mesh for the cells
        ctx[:contourplot] = XMakie.poly!(
            ctx[:scene],
            map(data -> make_mesh(data.g, data.f, data.e), ctx[:contourdata]);
            color = map(data -> vcat(data.f...), ctx[:contourdata]),
            colorrange = map(data -> data.c, ctx[:contourdata]),
            colormap = ctx[:colormap],
        )

        # draw the isolines via marching triangles
        if ctx[:elevation] ≈ 0
            XMakie.linesegments!(
                ctx[:scene],
                map(data -> first(marching_triangles(data.g, data.f, [], data.l; gridscale)), ctx[:contourdata]);
                color = :black,
                linewidth = ctx[:linewidth],
            )
        end

        XMakie.reset_limits!(ctx[:scene])
        add_scene!(ctx, makescene2d(ctx, :contourplot))
    end
    return reveal(ctx, TP)
end

# 2D vector
function vectorplot!(ctx, TP::Type{MakieType}, ::Type{Val{2}}, grid, func)
    XMakie = ctx[:Plotter]

    rc, rv = vectorsample(grid, func; gridscale = ctx[:gridscale], rasterpoints = ctx[:rasterpoints], offset = ctx[:offset])
    qc, qv = quiverdata(rc, rv; vscale = ctx[:vscale], vnormalize = ctx[:vnormalize], vconstant = ctx[:vconstant])

    set_plot_data!(ctx, :arrowdata, (qc = qc, qv = qv))

    if !haskey(ctx, :arrowplot)
        if !haskey(ctx, :scene)
            ctx[:scene] = XMakie.Axis(
                ctx[:figure];
                title = ctx[:title],
                aspect = XMakie.DataAspect(),
                scenekwargs(ctx)...,
            )
            add_scene!(ctx, ctx[:scene])
        end
        if isdefined(XMakie, :arrows2d!) # Makie >=0.23
            ctx[:arrowplot] = XMakie.arrows2d!(
                ctx[:scene],
                map(data -> data.qc[1, :], ctx[:arrowdata]),
                map(data -> data.qc[2, :], ctx[:arrowdata]),
                map(data -> data.qv[1, :], ctx[:arrowdata]),
                map(data -> data.qv[2, :], ctx[:arrowdata]);
                color = :black,
                shaftwidth = ctx[:linewidth],
                tipwidth = 3 * ctx[:linewidth],
            )
        else
            ctx[:arrowplot] = XMakie.arrows!(
                ctx[:scene],
                map(data -> data.qc[1, :], ctx[:arrowdata]),
                map(data -> data.qc[2, :], ctx[:arrowdata]),
                map(data -> data.qv[1, :], ctx[:arrowdata]),
                map(data -> data.qv[2, :], ctx[:arrowdata]);
                color = :black,
                linewidth = ctx[:linewidth],
            )
        end
        XMakie.reset_limits!(ctx[:scene])
    end
    return reveal(ctx, TP)
end

function streamplot!(ctx, TP::Type{MakieType}, ::Type{Val{2}}, grid, func)
    XMakie = ctx[:Plotter]

    rc, rv = vectorsample(
        grid, func; rasterpoints = 2 * ctx[:rasterpoints], offset = ctx[:offset], xlimits = ctx[:xlimits],
        ylimits = ctx[:ylimits], gridscale = ctx[:gridscale]
    )

    x = rc[1]
    y = rc[2]

    ix = linear_interpolation((x, y), rv[1, :, :])
    iy = linear_interpolation((x, y), rv[2, :, :])
    f(x, y) = Point2(ix(x, y), iy(x, y))

    xextent = x[end] - x[begin]
    yextent = y[end] - y[begin]

    maxextent = max(xextent, yextent)
    gridstep = maxextent / (2 * ctx[:rasterpoints])
    gridsize = (Int(ceil(xextent / gridstep)), Int(ceil(yextent / gridstep)), 2 * ctx[:rasterpoints])

    set_plot_data!(ctx, :streamdata, (xinterval = x[begin] .. x[end], yinterval = y[begin] .. y[end], f = f))

    if !haskey(ctx, :streamplot)
        if !haskey(ctx, :scene)
            ctx[:scene] = XMakie.Axis(
                ctx[:figure];
                title = ctx[:title],
                aspect = XMakie.DataAspect(),
                scenekwargs(ctx)...,
            )
            add_scene!(ctx, ctx[:scene])
        end
        ctx[:streamplot] = XMakie.streamplot!(
            ctx[:scene],
            map(data -> data.f, ctx[:streamdata]),
            map(data -> data.xinterval, ctx[:streamdata]),
            map(data -> data.yinterval, ctx[:streamdata]);
            linewidth = ctx[:linewidth],
            colormap = ctx[:colormap],
            gridsize = gridsize,
            arrow_size = 7.5,
            stepsize = 0.01 * maxextent
        )
        XMakie.reset_limits!(ctx[:scene])
    end
    return reveal(ctx, TP)
end

#######################################################################################
#######################################################################################
# 3D Grid

function xyzminmax(grid::ExtendableGrid, gridscale)
    coord = grid[Coordinates]
    ndim = size(coord, 1)
    xyzmin = zeros(ndim)
    xyzmax = ones(ndim)
    for idim in 1:ndim
        @views mn, mx = extrema(coord[idim, :])
        xyzmin[idim] = mn * gridscale
        xyzmax[idim] = mx * gridscale
    end
    return xyzmin, xyzmax
end

"""
   makeaxis3d(ctx)

Dispatch between LScene and new Axis3. Axis3 does not allow zoom, so we
support LScene in addition.
"""
function makeaxis3d(ctx)
    XMakie = ctx[:Plotter]
    if ctx[:scene3d] == :LScene
        # "Old" LScene with zoom-in functionality
        lim = ctx[:limits]
        lim = Observable(Rect3f(Vec3f(lim[1], lim[3], lim[5]), Vec3f(lim[2], lim[4], lim[6])))
        scene = XMakie.LScene(ctx[:figure], scenekw = (; limits = lim))
        return scene
    else
        # "New" Axis3 with prospective new stuff by Julius.
        return XMakie.Axis3(
            ctx[:figure];
            aspect = :data,
            viewmode = ctx[:viewmode],
            limits = ctx[:limits],
            elevation = ctx[:elev] * π / 180,
            azimuth = ctx[:azim] * π / 180,
            perspectiveness = ctx[:perspectiveness],
            title = map(data -> data.t, ctx[:data]),
            scenekwargs(ctx)...,
        )
    end
end

"""
       makescene3d(ctx)

Complete scene with title and status line showing interaction state.
This uses a gridlayout and its  protrusion capabilities.
"""
function makescene3d(ctx)
    XMakie = ctx[:Plotter]
    GL = XMakie.GridLayout(ctx[:figure]; default_rowgap = 0)
    if ctx[:scene3d] == "LScene"
        # LScene has no title, put the title into protrusion space on top  of the scene
        GL[1, 1, XMakie.Top()] = XMakie.Label(
            ctx[:figure],
            " $(map(data -> data.t, ctx[:data])) ";
            tellwidth = false,
            height = 30,
            fontsize = ctx[:fontsize],
        )
    end
    GL[1, 1] = ctx[:scene]
    # Horizontal or vertical colorbar
    if ctx[:show_colorbar]
        if haskey(ctx, :crange)
            if ctx[:colorbar] == :vertical
                GL[1, 2] = XMakie.Colorbar(
                    ctx[:figure];
                    colormap = ctx[:colormap],
                    colorrange = ctx[:crange],
                    ticks = map(d -> d.c, ctx[:data]),
                    tickformat = "{:.2e}",
                    width = 15,
                    ticklabelsize = 0.5 * ctx[:fontsize],
                )
            elseif ctx[:colorbar] == :horizontal
                GL[2, 1] = XMakie.Colorbar(
                    ctx[:figure];
                    colormap = ctx[:colormap],
                    colorrange = ctx[:crange],
                    ticks = map(d -> d.c, ctx[:data]),
                    tickformat = "{:.2e}",
                    height = 15,
                    ticklabelsize = 0.5 * ctx[:fontsize],
                    vertical = false,
                )
            end
        end
    end
    # Put the status label into protrusion space on the bottom of the scene
    GL[1, 1, XMakie.Bottom()] = XMakie.Label(
        ctx[:figure],
        ctx[:status];
        tellwidth = false,
        height = 30,
        fontsize = 0.5 * ctx[:fontsize],
    )
    return GL
end

const keyboardhelp = """
    Keyboard interactions:
          x: control xplane
          y: control yplane
          z: control zplane
          l: control isolevel
    up/down: fine control control value
pgup/pgdown: coarse control control value
          h: print this message
"""

function gridplot!(ctx, TP::Type{MakieType}, ::Type{Val{3}}, grid)
    function make_mesh(pts, fcs)
        if pkgversion(GeometryBasics) < v"0.5"
            return GeometryBasics.Mesh(meta(pts; normals = normals(pts, fcs)), fcs)
        else
            return GeometryBasics.Mesh(pts, fcs; normal = normals(pts, fcs))
        end
    end
    nregions = num_cellregions(grid)
    nbregions = num_bfaceregions(grid)

    XMakie = ctx[:Plotter]
    xyzmin, xyzmax = xyzminmax(grid, ctx[:gridscale])
    xyzstep = (xyzmax - xyzmin) / 100
    ctx[:limits] = (
        xyzmin[1], xyzmax[1],
        xyzmin[2], xyzmax[2],
        xyzmin[3], xyzmax[3],
    )

    function adjust_planes(xplane, yplane, zplane)
        ctx[:ixplane] = max(xyzmin[1], min(xyzmax[1], xplane))
        ctx[:iyplane] = max(xyzmin[2], min(xyzmax[2], yplane))
        return ctx[:izplane] = max(xyzmin[3], min(xyzmax[3], zplane))
    end

    adjust_planes(
        ctx[:xplanes][1],
        ctx[:yplanes][1],
        ctx[:zplanes][1]
    )

    if !haskey(ctx, :scene)
        ctx[:data] = Observable(
            (
                g = grid,
                x = ctx[:ixplane],
                y = ctx[:iyplane],
                z = ctx[:izplane],
                t = ctx[:title],
            )
        )

        ctx[:scene] = makeaxis3d(ctx)
        cmap = region_cmap(nregions)
        ctx[:cmap] = cmap
        bcmap = bregion_cmap(nbregions)
        ctx[:bcmap] = bcmap

        ############# Interior cuts
        # We draw a mesh for each color.
        if ctx[:interior]
            ctx[:celldata] = map(
                d -> extract_visible_cells3D(
                    d.g,
                    [d.x, d.y, d.z] / ctx[:gridscale];
                    cellcoloring = ctx[:cellcoloring],
                    gridscale = ctx[:gridscale],
                    primepoints = hcat(xyzmin, xyzmax),
                    Tp = Point3f,
                    Tf = GLTriangleFace,
                ),
                ctx[:data]
            )

            ctx[:cellmeshes] = map(d -> [make_mesh(d[1][i], d[2][i]) for i in 1:nregions], ctx[:celldata])

            for i in 1:nregions
                XMakie.mesh!(
                    ctx[:scene],
                    map(d -> d[i], ctx[:cellmeshes]);
                    color = cmap[i],
                    backlight = 1.0f0,
                )

                if ctx[:linewidth] > 0
                    XMakie.wireframe!(
                        ctx[:scene],
                        map(d -> d[i], ctx[:cellmeshes]);
                        color = :black,
                        linewidth = ctx[:linewidth],
                    )
                end
            end
        end

        ############# Visible boundary faces
        ctx[:facedata] = map(
            d -> extract_visible_bfaces3D(
                d.g,
                [d.x, d.y, d.z] / ctx[:gridscale];
                gridscale = ctx[:gridscale],
                primepoints = hcat(xyzmin, xyzmax),
                Tp = Point3f,
                Tf = GLTriangleFace,
            ),
            ctx[:data]
        )

        ctx[:facemeshes] = map(d -> [make_mesh(d[1][i], d[2][i]) for i in 1:nbregions], ctx[:facedata])

        for i in 1:nbregions
            XMakie.mesh!(
                ctx[:scene],
                map(d -> d[i], ctx[:facemeshes]);
                color = bcmap[i],
                backlight = 1.0f0,
            )
            if ctx[:linewidth] > 0
                XMakie.wireframe!(
                    ctx[:scene],
                    map(d -> d[i], ctx[:facemeshes]);
                    color = :black,
                    linewidth = ctx[:linewidth],
                )
            end
        end

        ############# Transparent outline

        if ctx[:outlinealpha] > 0.0
            ctx[:outlinedata] = map(
                d -> extract_visible_bfaces3D(
                    d.g,
                    xyzmax / ctx[:gridscale];
                    gridscale = ctx[:gridscale],
                    primepoints = hcat(xyzmin, xyzmax),
                    Tp = Point3f,
                    Tf = GLTriangleFace,
                ),
                ctx[:data]
            )
            ctx[:outlinemeshes] = map(
                d -> [make_mesh(d[1][i], d[2][i]) for i in 1:nbregions],
                ctx[:outlinedata]
            )

            for i in 1:nbregions
                XMakie.mesh!(
                    ctx[:scene],
                    map(d -> d[i], ctx[:outlinemeshes]);
                    color = (bcmap[i], ctx[:outlinealpha]),
                    transparency = true,
                    backlight = 1.0f0,
                )
            end
        end

        ##### Interaction
        scene_interaction(ctx[:scene].scene, XMakie, [:z, :y, :x, :q]) do delta, key
            if key == :x
                ctx[:ixplane] += delta * xyzstep[1]
                ctx[:status][] = @sprintf("x=%.3g", ctx[:ixplane])
            elseif key == :y
                ctx[:iyplane] += delta * xyzstep[2]
                ctx[:status][] = @sprintf("y=%.3g", ctx[:iyplane])
            elseif key == :z
                ctx[:izplane] += delta * xyzstep[3]
                ctx[:status][] = @sprintf("z=%.3g", ctx[:izplane])
            elseif key == :q
                ctx[:status][] = " "
            end

            adjust_planes(ctx[:ixplane], ctx[:iyplane], ctx[:izplane])

            ctx[:data][] = (
                g = grid,
                x = ctx[:ixplane],
                y = ctx[:iyplane],
                z = ctx[:izplane],
                t = ctx[:title],
            )
        end

        ctx[:status] = Observable(" ")

        XMakie.reset_limits!(ctx[:scene])
        add_scene!(ctx, makescene_grid(ctx))

    else
        ctx[:data][] = (
            g = grid,
            x = ctx[:ixplane],
            y = ctx[:iyplane],
            z = ctx[:izplane],
            t = ctx[:title],
        )
    end

    return reveal(ctx, TP)
end

# 3d function
function scalarplot!(ctx, TP::Type{MakieType}, ::Type{Val{3}}, grids, parentgrid, funcs)
    levels, crange, colorbarticks = isolevels(ctx, funcs)
    ctx[:crange] = crange
    ctx[:colorbarticks] = colorbarticks

    nan_replacement = 0.5 * (crange[1] + crange[2])
    make_mesh(pts, fcs) = Mesh(pts, fcs)

    function make_mesh(pts, fcs, vals, alpha)
        if length(fcs) > 0
            colors = XMakie.Makie.interpolated_getindex.((cmap,), vals, (crange,))
            if alpha < 1
                colors = [
                    RGBA(colors[i].r, colors[i].g, colors[i].b, Float32(alpha)) for
                        i in 1:length(colors)
                ]
            end
            if pkgversion(GeometryBasics) < v"0.5"
                return GeometryBasics.Mesh(meta(pts; color = colors, normals = normals(pts, fcs)), fcs)
            else
                return GeometryBasics.Mesh(pts, fcs; color = colors, normal = normals(pts, fcs))
            end
        else
            return GeometryBasics.Mesh(pts, fcs)
        end
    end

    nbregions = num_bfaceregions(parentgrid)

    XMakie = ctx[:Plotter]
    cmap = XMakie.to_colormap(ctx[:colormap])
    xyzmin, xyzmax = xyzminmax(parentgrid, ctx[:gridscale])
    xyzstep = (xyzmax - xyzmin) / 100

    fstep = (crange[2] - crange[1]) / 100
    if fstep ≈ 0
        fstep = 0.1
    end

    ctx[:limits] = (
        xyzmin[1], xyzmax[1],
        xyzmin[2], xyzmax[2],
        xyzmin[3], xyzmax[3],
    )

    ctx[:ixplanes] = collect(ctx[:xplanes]) * ctx[:gridscale]
    ctx[:iyplanes] = collect(ctx[:yplanes]) * ctx[:gridscale]
    ctx[:izplanes] = collect(ctx[:zplanes]) * ctx[:gridscale]

    x = ctx[:ixplanes]
    y = ctx[:iyplanes]
    z = ctx[:izplanes]

    ε = 1.0e-5 * (xyzmax .- xyzmin)

    ctx[:ixplanes] = isa(x, Number) ?
        collect(range(xyzmin[1] + ε[1], xyzmax[1] - ε[1]; length = ceil(x))) : x
    ctx[:iyplanes] = isa(y, Number) ?
        collect(range(xyzmin[2] + ε[2], xyzmax[2] - ε[2]; length = ceil(y))) : y
    ctx[:izplanes] = isa(z, Number) ?
        collect(range(xyzmin[3] + ε[3], xyzmax[3] - ε[3]; length = ceil(z))) : z

    ctx[:ixplanes] = max.(xyzmin[1], min.(xyzmax[1], ctx[:ixplanes]))
    ctx[:iyplanes] = max.(xyzmin[2], min.(xyzmax[2], ctx[:iyplanes]))
    ctx[:izplanes] = max.(xyzmin[3], min.(xyzmax[3], ctx[:izplanes]))

    ctx[:levels] = levels

    if !haskey(ctx, :scene)
        ctx[:data] = Observable(
            (
                g = grids,
                p = parentgrid,
                f = funcs,
                x = ctx[:ixplanes],
                y = ctx[:iyplanes],
                z = ctx[:izplanes],
                l = ctx[:levels],
                c = ctx[:colorbarticks],
                t = ctx[:title],
            )
        )

        ctx[:scene] = makeaxis3d(ctx)

        #### Transparent outline
        if ctx[:outlinealpha] > 0.0
            ctx[:outlinedata] = map(
                d -> extract_visible_bfaces3D(
                    d.p,
                    xyzmax / ctx[:gridscale];
                    gridscale = ctx[:gridscale],
                    primepoints = hcat(xyzmin, xyzmax),
                    Tp = Point3f,
                    Tf = GLTriangleFace,
                ),
                ctx[:data]
            )
            ctx[:facemeshes] = map(
                d -> [make_mesh(d[1][i], d[2][i]) for i in 1:nbregions],
                ctx[:outlinedata]
            )
            bcmap = bregion_cmap(nbregions)
            for i in 1:nbregions
                XMakie.mesh!(
                    ctx[:scene],
                    map(d -> d[i], ctx[:facemeshes]);
                    color = (bcmap[i], ctx[:outlinealpha]),
                    transparency = true,
                    backlight = 1.0f0,
                )
            end
        end

        make_planes = d -> make_mesh(
            marching_tetrahedra(
                d.g,
                d.f,
                makeplanes(xyzmin, xyzmax, d.x, d.y, d.z),
                [];
                gridscale = ctx[:gridscale],
                primepoints = hcat(xyzmin, xyzmax),
                primevalues = crange,
                tol = ctx[:tetxplane_tol],
                Tp = Point3f,
                Tf = GLTriangleFace,
                Tv = Float32,
            )...,
            ctx[:planealpha]
        )

        make_levels = d -> make_mesh(
            marching_tetrahedra(
                d.g,
                d.f,
                [],
                d.l;
                gridscale = ctx[:gridscale],
                primepoints = hcat(xyzmin, xyzmax),
                primevalues = crange,
                tol = ctx[:tetxplane_tol],
                Tp = Point3f,
                Tf = GLTriangleFace,
                Tv = Float32,
            )...,
            ctx[:levelalpha]
        )

        #### Plane sections and isosurfaces
        ctx[:planesections] = XMakie.mesh!(
            ctx[:scene],
            map(make_planes, ctx[:data]);
            backlight = 1.0f0,
            transparency = ctx[:planealpha] < 1.0,
        )

        ctx[:isosurfaces] = XMakie.mesh!(
            ctx[:scene],
            map(make_levels, ctx[:data]);
            backlight = 1.0f0,
            transparency = ctx[:levelalpha] < 1.0,
        )

        #### Interactions
        scene_interaction(ctx[:scene].scene, XMakie, [:z, :y, :x, :l, :q]) do delta, key
            if key == :x
                ctx[:ixplanes] .+= delta * xyzstep[1]
                ctx[:status][] = "x=[" * mapreduce(x -> @sprintf("%.3g,", x), *, ctx[:ixplanes][1]) * "]"
            elseif key == :y
                ctx[:iyplanes] .+= delta * xyzstep[2]
                ctx[:status][] = "y=[" * mapreduce(y -> @sprintf("%.3g,", y), *, ctx[:iyplanes][1]) * "]"
            elseif key == :z
                ctx[:izplanes] .+= delta * xyzstep[3]
                ctx[:status][] = "z=[" * mapreduce(z -> @sprintf("%.3g,", z), *, ctx[:izplanes][1]) * "]"
            elseif key == :l
                ctx[:levels] .+= delta * fstep
                ctx[:status][] = "l=[" * mapreduce(l -> @sprintf("%.3g,", l), *, ctx[:levels]) * "]"
            elseif key == :q
                ctx[:status][] = " "
            end

            ctx[:ixplanes] = max.(xyzmin[1], min.(xyzmax[1], ctx[:ixplanes]))
            ctx[:iyplanes] = max.(xyzmin[2], min.(xyzmax[2], ctx[:iyplanes]))
            ctx[:izplanes] = max.(xyzmin[3], min.(xyzmax[3], ctx[:izplanes]))

            ctx[:data][] = (
                g = grids,
                p = parentgrid,
                f = funcs,
                x = ctx[:ixplanes],
                y = ctx[:iyplanes],
                z = ctx[:izplanes],
                l = ctx[:levels],
                t = ctx[:title],
            )
        end
        ctx[:status] = Observable(" ")
        add_scene!(ctx, makescene3d(ctx))
    else
        ctx[:data][] = (
            g = grids,
            p = parentgrid,
            f = funcs,
            x = ctx[:ixplanes],
            y = ctx[:iyplanes],
            z = ctx[:izplanes],
            l = ctx[:levels],
            t = ctx[:title],
        )
    end
    return reveal(ctx, TP)
end

# TODO: allow aspect scaling
# if ctx[:aspect]>1.0
#     XMakie.scale!(ctx[:scene],ctx[:aspect],1.0)
# else
#     XMakie.scale!(ctx[:scene],1.0,1.0/ctx[:aspect])
# end

# TODO: use distinguishable colors
# http://juliagraphics.github.io/Colors.jl/stable/colormapsandcolorscales/#Generating-distinguishable-colors-1

# TODO: a priori angles aka pyplot3D
# rect = ctx[:scene]
# azim=ctx[:azim]
# elev=ctx[:elev]
# arr = normalize([cosd(azim/2), 0, sind(azim/2), -sind(azim/2)])
# XMakie.rotate!(rect, XMakie.Quaternionf0(arr...))

# 3 replies
# Julius Krumbiegel  5 hours ago
# try lines(x, y, axis = (limits = lims,)), the other keyword arguments go to the plot
# Julius Krumbiegel  5 hours ago
# although I think it should be targetlimits because the limits are computed from them usually, considering that there might be aspect constraints that should be met or linked axes (edited)
# Christophe Meyer  4 hours ago
# Thanks!  lines(x, y, axis = (targetlimits = lims,))  indeed makes the limits update.^
# I found that autolimits!(axis) gave good results, even better than me manually computing limits!

function customplot!(ctx, TP::Type{MakieType}, func)
    XMakie = ctx[:Plotter]
    if !haskey(ctx, :scene)
        ctx[:scene] = XMakie.Axis(
            ctx[:figure];
            title = ctx[:title],
            aspect = XMakie.DataAspect(),
            scenekwargs(ctx)...,
        )
        add_scene!(ctx, ctx[:scene])
    end
    func(ctx[:scene])
    return reveal(ctx, TP)
end

#
# Code moved from https://github.com/JuliaGeometry/Triangulate.jl/blob/2aa1196625832b622eb808b165cab8c7a098c6c5/src/plot.jl
#
function plot_triangulateio!(
        ctx,
        TP::Type{MakieType},
        triangulateio;
        voronoi = nothing,
        circumcircles = false,
        kwargs...
    )
    XMakie = ctx[:Plotter]

    function frgb(Plotter, i, max; pastel = false)
        x = Float64(i - 1) / Float64(max)
        if (x < 0.5)
            r = 1.0 - 2.0 * x
            g = 2.0 * x
            b = 0.0
        else
            r = 0.0
            g = 2.0 - 2.0 * x
            b = 2.0 * x - 1.0
        end
        if pastel
            r = 0.5 + 0.5 * r
            g = 0.5 + 0.5 * g
            b = 0.5 + 0.5 * b
        end
        return XMakie.RGBA(r, g, b)
    end


    x = view(triangulateio.pointlist, 1, :)
    y = view(triangulateio.pointlist, 2, :)

    xminmax = extrema(x)
    yminmax = extrema(y)

    dx = xminmax[2] - xminmax[1]
    dy = yminmax[2] - yminmax[1]

    if !haskey(ctx, :scene)
        if ctx[:aspect] ≈ 1.0
            aspect = XMakie.DataAspect()
        else
            autolimitaspect = ctx[:aspect]
        end
        ctx[:scene] = XMakie.Axis(
            ctx[:figure];
            title = ctx[:title],
            aspect,
            scenekwargs(ctx)...,
        )
    end
    axis = ctx[:scene]


    if size(triangulateio.pointlist, 2) > 0
        XMakie.ylims!(axis, (yminmax[1] - dy / 10, yminmax[2] + dy / 10))
        XMakie.xlims!(axis, (xminmax[1] - dx / 10, xminmax[2] + dx / 10))

        points = reshape(reinterpret(XMakie.Point{2, Cdouble}, triangulateio.pointlist), size(triangulateio.pointlist, 2))
        if size(triangulateio.trianglelist, 2) > 0
            if size(triangulateio.triangleattributelist, 2) > 0
                attr = triangulateio.triangleattributelist[1, :]
                maxattr = maximum(attr)
                for i in 1:size(triangulateio.trianglelist, 2)
                    XMakie.poly!(
                        axis,
                        view(points, view(triangulateio.trianglelist, :, i));
                        strokecolor = :black,
                        strokewidth = 1,
                        color = frgb(XMakie, attr[i], maxattr; pastel = true)
                    )
                end

            else
                for i in 1:size(triangulateio.trianglelist, 2)
                    XMakie.poly!(
                        axis,
                        view(points, view(triangulateio.trianglelist, :, i));
                        strokecolor = :black,
                        strokewidth = 1,
                        color = :white,
                        alpha = 0.75
                    )
                end
            end
        end

        if circumcircles
            col = XMakie.RGB(0.4, 0.05, 0.4)

            t = 0:(0.025 * π):(2π)
            circle(x, y, r) = XMakie.lines!(
                axis,
                x .+ r .* cos.(t),
                y .+ r .* sin.(t);
                color = col,
                linewidth = 0.3
            )
            cc = zeros(2)
            for itri in 1:size(triangulateio.trianglelist, 2)
                a = view(triangulateio.pointlist, :, triangulateio.trianglelist[1, itri])
                b = view(triangulateio.pointlist, :, triangulateio.trianglelist[2, itri])
                c = view(triangulateio.pointlist, :, triangulateio.trianglelist[3, itri])
                tricircumcenter!(cc, a, b, c)
                r = sqrt((cc[1] - a[1])^2 + (cc[2] - a[2])^2)
                circle(cc[1], cc[2], r)
                XMakie.scatter!(axis, [cc[1]], [cc[2]]; markersize = 5, color = col)
            end
        end
        if size(triangulateio.segmentlist, 2) > 0
            markermax = maximum(triangulateio.segmentmarkerlist)
            # see https://gist.github.com/gizmaa/7214002
            xx = zeros(2)
            yy = zeros(2)
            for i in 1:size(triangulateio.segmentlist, 2)
                xx[1] = triangulateio.pointlist[1, triangulateio.segmentlist[1, i]]
                yy[1] = triangulateio.pointlist[2, triangulateio.segmentlist[1, i]]
                xx[2] = triangulateio.pointlist[1, triangulateio.segmentlist[2, i]]
                yy[2] = triangulateio.pointlist[2, triangulateio.segmentlist[2, i]]
                XMakie.lines!(
                    axis,
                    xx,
                    yy;
                    color = frgb(XMakie, triangulateio.segmentmarkerlist[i], markermax),
                    linewidth = 3
                )
            end
        end

        if size(triangulateio.trianglelist, 2) == 0 && size(triangulateio.regionlist, 2) > 0
            markermax = maximum(triangulateio.regionlist)
            xx = triangulateio.regionlist[1, :]
            yy = triangulateio.regionlist[2, :]
            reg = triangulateio.regionlist[3, :]
            color = [frgb(XMakie, reg[i], markermax) for i in 1:size(triangulateio.regionlist, 2)]
            XMakie.scatter!(
                axis,
                xx,
                yy;
                markersize = 10,
                color
            )
        end

        if size(triangulateio.trianglelist, 2) == 0 && size(triangulateio.holelist, 2) > 0
            xx = triangulateio.holelist[1, :]
            yy = triangulateio.holelist[2, :]
            XMakie.scatter!(axis, xx, yy; markersize = 10, color = :brown)
        end

        if voronoi != nothing && size(voronoi.edgelist, 2) > 0
            bcx = sum(x) / size(triangulateio.pointlist, 2)
            bcy = sum(y) / size(triangulateio.pointlist, 2)
            wx = maximum(abs.(x .- bcx))
            wy = maximum(abs.(y .- bcy))
            ww = max(wx, wy)

            x = voronoi.pointlist[1, :]
            y = voronoi.pointlist[2, :]
            XMakie.scatter!(axis, x, y; markersize = 10, color = :green)
            for i in 1:size(voronoi.edgelist, 2)
                i1 = voronoi.edgelist[1, i]
                i2 = voronoi.edgelist[2, i]
                if i1 > 0 && i2 > 0
                    XMakie.lines!(
                        axis,
                        [voronoi.pointlist[1, i1], voronoi.pointlist[1, i2]],
                        [voronoi.pointlist[2, i1], voronoi.pointlist[2, i2]];
                        color = :green
                    )
                else
                    x0 = voronoi.pointlist[1, i1]
                    y0 = voronoi.pointlist[2, i1]
                    xn = voronoi.normlist[1, i]
                    yn = voronoi.normlist[2, i]
                    normscale = 1.0e10
                    if x0 + normscale * xn > bcx + ww
                        normscale = min(normscale, abs((bcx + ww - x0) / xn))
                    end
                    if x0 + normscale * xn < bcx - ww
                        normscale = min(normscale, abs((bcx - ww - x0) / xn))
                    end
                    if y0 + normscale * yn > bcy + ww
                        normscale = min(normscale, abs((bcy + ww - y0) / yn))
                    end
                    if y0 + normscale * yn < bcy - ww
                        normscale = min(normscale, abs((bcy - ww - y0) / yn))
                    end
                    XMakie.lines!(
                        axis,
                        [x0, x0 + normscale * xn],
                        [y0, y0 + normscale * yn];
                        color = :green
                    )
                end
            end
        end
        XMakie.scatter!(axis, points; markersize = 7.5, color = :black)
        GL = XMakie.GridLayout(ctx[:figure])
        GL[1, 1] = ctx[:scene]
        add_scene!(ctx, GL)
    end
    return reveal(ctx, TP)
end


end
