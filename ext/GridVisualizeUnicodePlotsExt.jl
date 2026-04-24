"""
    module GridVisualizeUnicodePlotsExt
    
Extension module for UnicodePlots.jl
"""
module GridVisualizeUnicodePlotsExt

import GridVisualize: initialize!, gridplot!, scalarplot!, vectorplot!, bregion_cmap, region_cmap, reveal, streamplot!
using GridVisualize: UnicodePlotsType, GridVisualizer, SubVisualizer, vectorsample, quiverdata
using UnicodePlots: UnicodePlots
using ExtendableGrids: Coordinates, simplexgrid, ON_CELLS, ON_FACES, ON_EDGES, CellNodes, FaceNodes, BFaceNodes, CellGeometries, CellRegions, BFaceRegions, num_cells, num_nodes, local_celledgenodes, num_bfaceregions, num_cellregions, num_targets, interpolate!
using Colors: Colors, RGB, RGBA, red, green, blue
using ColorSchemes: colorschemes, color

initialize!(p, ::Type{UnicodePlotsType}) = nothing


function reveal(p::GridVisualizer, ::Type{UnicodePlotsType})
    layout = p.context[:layout]
    subplots = @views permutedims(p.subplots)[:]

    if layout == (1, 1)
        display(subplots[1][:figure])
    else
        if :Term ∉ Symbol.(Base.loaded_modules_array())
            @warn "A GridVisualizer with multiple UnicodePlots requires 'Term.jl' to be loaded: add Term.jl to your environment."
        else
            figures = [subplot[:figure] for subplot in subplots if haskey(subplot, :figure)]
            grid_plot = UnicodePlots.gridplot(figures, layout = p.context[:layout], show_placeholder = true)
            display(grid_plot)
        end
    end
    return nothing
end


function reveal(ctx::SubVisualizer, TP::Type{UnicodePlotsType})
    if ctx[:show] || ctx[:reveal]
        return reveal(ctx[:GridVisualizer], TP)
    end
    return nothing
end


function region_legend!(plt, title, y0, colors)
    # legend by annotate!
    UnicodePlots.label!(plt, :r, y0, title)
    for r in 1:length(colors)
        red, green, blue = UInt32.(colors[r])
        uint_color = (red << 16) | (green << 8) | blue
        y0 += 1
        UnicodePlots.label!(plt, :r, y0, "   " * string(r), uint_color)
    end
    return y0
end

gridplot!(ctx, TP::Type{UnicodePlotsType}, ::Type{Val{3}}, grid) = @warn "3D gridplots are not implemented for the UnicodePlots backend"

function gridplot!(ctx, TP::Type{UnicodePlotsType}, ::Type{Val{2}}, grid)
    UnicodePlots = ctx[:Plotter]

    # find bounding box
    coords = grid[Coordinates]
    xlimits = ctx[:xlimits]
    ylimits = ctx[:ylimits]
    if xlimits[1] < xlimits[2]
        ex = xlimits
    else
        ex = extrema(view(coords, 1, :))
    end
    if ylimits[1] < ylimits[2]
        ey = ylimits
    else
        ey = extrema(view(coords, 2, :))
    end

    # line color for interior edges
    edge_color = UnicodePlots.ansi_color(:normal)

    # determine resolution (divided by 10, to reduce pixel count in the terminal)
    layout = ctx[:layout]
    resolution = ctx[:size] ./ 12 ./ (layout[2], layout[1])
    aspect = ctx[:aspect]

    if (true) # auto scale feature, do we want this?
        wx = ex[2] - ex[1]
        wy = ey[2] - ey[1]
        rescale = wx / wy * (resolution[1] / (2 * resolution[2]))
        if rescale > 1
            resolution = (resolution[1] * aspect, Int(ceil(resolution[2] / rescale)))
        else
            resolution = (Int(ceil(resolution[1] * aspect / rescale)), resolution[2])
        end
    end

    # we need an integer resolution
    resolution = @. Int(round(resolution))

    # ensure that legend fits
    ncellregions = num_cellregions(grid)
    nbregions = num_bfaceregions(grid)
    resolution = (resolution[1], max(resolution[2], 5 + ncellregions + nbregions))

    # create UnicodePlots.Canvas
    CanvasType = UnicodePlots.BrailleCanvas # should this be a changeable parameter ?
    canvas = CanvasType(
        resolution[2], resolution[1],            # number of rows and columns (characters)
        origin_y = ey[1], origin_x = ex[1],      # position in virtual space
        height = (ey[2] - ey[1]), width = ex[2] - ex[1]; blend = false
    )

    linewidth = ctx[:linewidth]
    if linewidth > 0
        ## plot all edges in the grid
        plot_based = ctx[:cellwise] ? ON_CELLS : ON_FACES
        if plot_based in [ON_FACES, ON_EDGES]
            # plot all edges via FaceNodes
            facenodes = grid[FaceNodes]
            nfaces = size(facenodes, 2)
            for j in 1:nfaces
                UnicodePlots.lines!(
                    canvas,
                    coords[1, facenodes[1, j]], coords[2, facenodes[1, j]], # from
                    coords[1, facenodes[2, j]], coords[2, facenodes[2, j]]; # to
                    color = edge_color
                )
            end
        elseif plot_based == ON_CELLS
            # plot all edges via CellNodes and local_celledgenodes
            cellnodes = grid[CellNodes]
            cellgeoms = grid[CellGeometries]
            ncells = num_cells(grid)
            for j in 1:ncells
                cen = local_celledgenodes(cellgeoms[j])
                for k in 1:size(cen, 2)
                    UnicodePlots.lines!(
                        canvas,
                        coords[1, cellnodes[cen[1, k], j]], coords[2, cellnodes[cen[1, k], j]],
                        coords[1, cellnodes[cen[2, k], j]], coords[2, cellnodes[cen[2, k], j]];
                        color = color
                    )
                end
            end
        end
    end

    # color cell midpoints with cell regions color
    cellregions = grid[CellRegions]
    cmap = region_cmap(max(2, ncellregions))
    ctx[:cmap] = cmap
    cell_colors = [
        (
                Int(round(cmap[i].r * 255)),
                Int(round(cmap[i].g * 255)),
                Int(round(cmap[i].b * 255)),
            ) for i in 1:ncellregions
    ]
    cellnodes = grid[CellNodes]
    cellgeoms = grid[CellGeometries]
    ncells = num_cells(grid)
    midpoint = [0.0, 0.0]
    markersize = ctx[:markersize]
    for j in 1:ncells
        fill!(midpoint, 0.0)
        nvertices = num_targets(cellnodes, j)
        for k in 1:nvertices
            midpoint .+= coords[:, cellnodes[k, j]]
        end
        midpoint ./= nvertices
        r = cellregions[j]
        if markersize > 0
            if markersize < 4
                UnicodePlots.points!(
                    canvas,
                    midpoint[1], midpoint[2];
                    color = cell_colors[r]
                )
            else
                if markersize < 6
                    character = "•"
                elseif markersize < 8
                    character = "●"
                else
                    character = "⬤"
                end
                UnicodePlots.annotate!(
                    canvas,
                    midpoint[1], midpoint[2],
                    character,
                    UnicodePlots.ansi_color(cell_colors[r]),
                    false
                )
            end
        end

    end

    # plot boundary faces with bregion_cmap colors
    bcmap = bregion_cmap(nbregions)
    ctx[:bcmap] = bcmap
    bcolors = [
        (
                Int(round(bcmap[i].r * 255)),
                Int(round(bcmap[i].g * 255)),
                Int(round(bcmap[i].b * 255)),
            ) for i in 1:nbregions
    ]
    bfacenodes = grid[BFaceNodes]
    bfaceregions = grid[BFaceRegions]
    nbfaces = size(bfacenodes, 2)
    for j in 1:nbfaces
        UnicodePlots.lines!(
            canvas,
            coords[1, bfacenodes[1, j]], coords[2, bfacenodes[1, j]],
            coords[1, bfacenodes[2, j]], coords[2, bfacenodes[2, j]];
            color = bcolors[bfaceregions[j]]
        )
    end


    plt = UnicodePlots.Plot(canvas; title = ctx[:title], border = ctx[:border])

    y0 = 0
    if markersize > 0
        y0 = region_legend!(plt, " cell", 1, [])
        y0 = region_legend!(plt, "regions", 2, cell_colors)
    end
    region_legend!(plt, " bface", y0 + 2, [])
    region_legend!(plt, "regions", y0 + 3, bcolors)

    # corner coordinates
    UnicodePlots.label!(plt, :b, ctx[:xlabel])
    UnicodePlots.label!(plt, :bl, UnicodePlots.nice_repr(ex[1], plt), UnicodePlots.ansi_color(UnicodePlots.BORDER_COLOR[]))
    UnicodePlots.label!(plt, :br, UnicodePlots.nice_repr(ex[2], plt), UnicodePlots.ansi_color(UnicodePlots.BORDER_COLOR[]))
    UnicodePlots.label!(plt, :l, 1, UnicodePlots.nice_repr(ey[2], plt), UnicodePlots.ansi_color(UnicodePlots.BORDER_COLOR[]))
    UnicodePlots.label!(plt, :l, round(Int, (resolution[2] + 1) / 2), ctx[:ylabel])
    UnicodePlots.label!(plt, :l, resolution[2], UnicodePlots.nice_repr(ey[1], plt), UnicodePlots.ansi_color(UnicodePlots.BORDER_COLOR[]))

    # plot
    ctx[:figure] = plt
    return reveal(ctx, TP)
end


function gridplot!(ctx, TP::Type{UnicodePlotsType}, ::Type{Val{1}}, grid)
    UnicodePlots = ctx[:Plotter]

    # find bounding box
    xlimits = ctx[:xlimits]
    coords = grid[Coordinates]
    if xlimits[1] < xlimits[2]
        ex = xlimits
    else
        ex = extrema(view(coords, 1, :))
    end

    # determine resolution (divided by 5, to reduce pixel count in the terminal)
    ncellregions = num_cellregions(grid)
    nbregions = num_bfaceregions(grid)
    layout = ctx[:layout]
    resolution = (Int(round(ctx[:size][1] / 6 / layout[2])), max(7, 5 + ncellregions + nbregions))

    # create UnicodePlots.Canvas
    CanvasType = UnicodePlots.BrailleCanvas # should this be a changeable parameter ?
    canvas = CanvasType(
        resolution[2], resolution[1],               # number of rows and columns (characters)
        origin_y = 0, origin_x = ex[1],             # position in virtual space
        height = 1, width = ex[2] - ex[1]; blend = false
    )

    # plot all edges in the grid
    cellregions = grid[CellRegions]
    cmap = region_cmap(max(2, ncellregions))
    ctx[:cmap] = cmap
    colors = [
        (
                Int(round(cmap[i].r * 255)),
                Int(round(cmap[i].g * 255)),
                Int(round(cmap[i].b * 255)),
            ) for i in 1:ncellregions
    ]
    cellnodes = grid[CellNodes]
    cellgeoms = grid[CellGeometries]
    ncells = num_cells(grid)
    nnodes = num_nodes(grid)
    text_color = UnicodePlots.ansi_color(:normal)
    if nnodes < resolution[1] / 2
        for j in 1:nnodes
            UnicodePlots.annotate!(canvas, coords[1, j], 0.5, "•", text_color, false)
        end
    end
    for j in 1:ncells
        cen = local_celledgenodes(cellgeoms[j])
        r = cellregions[j]
        for k in 1:size(cen, 2)
            UnicodePlots.lines!(
                canvas,
                coords[1, cellnodes[cen[1, k], j]], 0.5,
                coords[1, cellnodes[cen[2, k], j]], 0.5;
                color = colors[r]
            )
        end
    end

    # plot boundary nodes with bregion_cmap colors
    bcmap = bregion_cmap(nbregions)
    ctx[:bcmap] = bcmap
    bcolors = [
        (
                Int(round(bcmap[i].r * 255)),
                Int(round(bcmap[i].g * 255)),
                Int(round(bcmap[i].b * 255)),
            ) for i in 1:nbregions
    ]
    bfacenodes = grid[BFaceNodes]
    bfaceregions = grid[BFaceRegions]
    nbfaces = size(bfacenodes, 2)
    for j in 1:nbfaces
        red, green, blue = UInt32.(bcolors[bfaceregions[j]])
        uint_color = (red << 16) | (green << 8) | blue
        UnicodePlots.annotate!(canvas, coords[1, bfacenodes[1, j]], 0.5, "•", uint_color, false)
    end

    plt = UnicodePlots.Plot(canvas; title = ctx[:title], border = ctx[:border])

    y0 = region_legend!(plt, " cell", 1, [])
    y0 = region_legend!(plt, "regions", 2, colors)
    region_legend!(plt, " bface", y0 + 2, [])
    region_legend!(plt, "regions", y0 + 3, bcolors)


    # corner coordinates
    UnicodePlots.label!(plt, :bl, string(Float16(ex[1])), UnicodePlots.ansi_color(UnicodePlots.BORDER_COLOR[]))
    UnicodePlots.label!(plt, :b, ctx[:xlabel])
    UnicodePlots.label!(plt, :br, string(Float16(ex[2])), UnicodePlots.ansi_color(UnicodePlots.BORDER_COLOR[]))

    # plot
    ctx[:figure] = plt

    return reveal(ctx, TP)
end


function scalarplot!(
        ctx,
        TP::Type{UnicodePlotsType},
        ::Type{Val{1}},
        grids,
        parentgrid,
        funcs
    )

    nfuncs = length(funcs)
    layout = ctx[:layout]
    resolution = @. Int(round(ctx[:size] ./ (12, 6) ./ (layout[2], 4 * layout[1]))) # reduce pixel count in the terminal (size is then compatible to other plots)
    @info resolution

    ylim = ctx[:limits]
    if ylim[1] > ylim[2]
        # try to find limits automatically
        ylim = (minimum([minimum(func) for func in funcs]), maximum([maximum(func) for func in funcs]))
    end

    if ctx[:clear] || !haskey(ctx, :figure)
        plt = nothing
    else
        plt = ctx[:figure]
    end

    xlim = ctx[:xlimits]
    if xlim[1] > xlim[2]
        # invalid, try to find the optimal range
        coord_min = min(ctx[:xlimits][1], minimum.([grid[Coordinates] for grid in grids])...)
        coord_max = max(ctx[:xlimits][2], maximum.([grid[Coordinates] for grid in grids])...)
        xlim = (coord_min, coord_max)
    end

    xscale = ctx[:xscale]
    xscale == :log && (xscale = :log10)
    xscale == :symlog && (xscale = x -> sign(x) * (log10(1 + abs(x))))

    yscale = ctx[:yscale]
    yscale == :log && (yscale = :log10)
    yscale == :symlog && (yscale = x -> sign(x) * (log10(1 + abs(x))))

    if typeof(ctx[:color]) <: String || typeof(ctx[:color]) <: Symbol
        color = UnicodePlots.ansi_color(Symbol(ctx[:color]))
    elseif typeof(ctx[:color]) <: RGB
        color = (
            Int(round(ctx[:color].r * 255)),
            Int(round(ctx[:color].g * 255)),
            Int(round(ctx[:color].b * 255)),
        )
    else
        color = ctx[:color]
    end

    for ifunc in 1:nfuncs
        func = funcs[ifunc]
        grid = grids[ifunc]
        coord = grid[Coordinates] * ctx[:gridscale]
        name = name = isnothing(ctx[:label]) ? "" : ctx[:label]

        if isnothing(plt)
            plt = UnicodePlots.lineplot(
                coord[1, :],
                func;
                xlim,
                ylim,
                xscale,
                yscale,
                xlabel = String(ctx[:xlabel]),
                ylabel = ctx[:ylabel],
                name,
                height = resolution[2],
                width = resolution[1],
                title = ctx[:title],
                border = ctx[:border],
                color
            )
        else
            UnicodePlots.lineplot!(
                plt,
                coord[1, :],
                func;
                name,
                color
            )
        end
    end

    ctx[:figure] = plt

    return reveal(ctx, TP)
end


function scalarplot!(
        ctx,
        TP::Type{UnicodePlotsType},
        ::Type{Val{2}},
        grids,
        parentgrid,
        funcs
    )

    func = funcs[1]
    layout = ctx[:layout]
    resolution = ctx[:size] ./ (12, 6) ./ (layout[2], layout[1]) # reduce pixel count in the terminal
    ylim = ctx[:limits]
    colormap = ctx[:colormap]

    if ylim[1] > ylim[2]
        # try to find limits automatically
        ylim = (minimum([minimum(func) for func in funcs]), maximum([maximum(func) for func in funcs]))
    end

    coords = grids[1][Coordinates]
    ex = extrema(view(coords, 1, :))
    ey = extrema(view(coords, 2, :))

    aspect = ctx[:aspect]
    if (true) # auto scale feature, do we want this?
        wx = ex[2] - ex[1]
        wy = ey[2] - ey[1]
        rescale = wx / wy * (resolution[1] / (resolution[2]))
        if rescale > 1
            resolution = (resolution[1] * aspect, Int(ceil(resolution[2] / rescale)))
        else
            resolution = (Int(ceil(resolution[1] * aspect / rescale)), resolution[2])
        end
    end

    # we need an integer resolution
    resolution = @. Int(round(resolution))

    X = LinRange(ex[1], ex[2], resolution[1])
    Y = LinRange(ey[1], ey[2], resolution[2])
    xgrid_plot = simplexgrid(X, Y)

    # interpolate data onto plot_grid
    I = zeros(Float64, num_nodes(xgrid_plot))
    interpolate!(I, xgrid_plot, func, grids[1]; eps = 1.0e-14, not_in_domain_value = NaN, trybrute = true)

    plt = UnicodePlots.heatmap(
        reshape(I, (resolution[1], resolution[2]))',
        xlabel = ctx[:xlabel],
        ylabel = ctx[:ylabel],
        xfact = (ex[2] - ex[1]) / (resolution[1] - 1),
        yfact = (ey[2] - ey[1]) / (resolution[2] - 1),
        xoffset = ex[1],
        yoffset = ey[1],
        title = ctx[:title],
        colormap = colormap,
        height = resolution[2],
        width = resolution[1],
        xscale = ctx[:xscale],
        yscale = ctx[:yscale],
        border = ctx[:border]
    )

    ctx[:figure] = plt

    return reveal(ctx, TP)
end

scalarplot!(ctx, TP::Type{UnicodePlotsType}, ::Type{Val{3}}, grids, parentgrid, funcs) = @warn "3D scalarplot is not implemented for the UnicodePlots backend"


# unicode arrows for vector plot
arrows_verythin = ['↙', '↓', '↘', '→', '↗', '↑', '↖', '←']
arrows_thin = ['🡯', '🡫', '🡮', '🡪', '🡭', '🡩', '🡬', '🡨']
arrows_medium = ['🡷', '🡳', '🡶', '🡲', '🡵', '🡱', '🡴', '🡰']
arrows_thick = ['🡿', '🡻', '🡾', '🡺', '🡽', '🡹', '🡼', '🡸']
arrows_verythick = ['🢇', '🢃', '🢆', '🢂', '🢅', '🢁', '🢄', '🢀']

# helper function that selects the right arrow for a given vector direction and norm
function select_arrow(angle, norm, scale)
    if norm * scale < 1.0e-2
        return '•' # use a dot for very small vectors
    end
    if angle > -7 * π / 8 && angle <= -5 * π / 8
        a = 1
    elseif angle > -5 * π / 8 && angle <= -3 * π / 8
        a = 2
    elseif angle > -3 * π / 8 && angle <= -π / 8
        a = 3
    elseif angle > -π / 8 && angle <= π / 8
        a = 4
    elseif angle > π / 8 && angle <= 3 * π / 8
        a = 5
    elseif angle > 3 * π / 8 && angle <= 5 * π / 8
        a = 6
    elseif angle > 5 * π / 8 && angle <= 7 * π / 8
        a = 7
    else
        a = 8
    end
    if norm * scale <= 0.2
        return arrows_verythin[a]
    elseif norm * scale <= 0.4
        return arrows_thin[a]
    elseif norm * scale <= 0.6
        return arrows_medium[a]
    elseif norm * scale <= 0.8
        return arrows_thick[a]
    else
        return arrows_verythick[a]
    end
end


function vectorplot!(ctx, TP::Type{UnicodePlotsType}, ::Type{Val{2}}, grid, func)

    layout = ctx[:layout]
    resolution = ctx[:size] ./ 12 ./ (layout[2], layout[1]) # reduce pixel count in the terminal

    # find bounding box
    coords = grid[Coordinates]
    xlimits = ctx[:xlimits]
    ylimits = ctx[:ylimits]
    if xlimits[1] < xlimits[2]
        ex = xlimits
    else
        ex = extrema(view(coords, 1, :))
    end
    if ylimits[1] < ylimits[2]
        ey = ylimits
    else
        ey = extrema(view(coords, 2, :))
    end
    aspect = ctx[:aspect]

    if (true) # auto scale feature, do we want this?
        wx = ex[2] - ex[1]
        wy = ey[2] - ey[1]
        rescale = wx / wy * (resolution[1] / (2 * resolution[2]))
        if rescale > 1
            resolution = (resolution[1] * aspect, Int(ceil(resolution[2] / rescale)))
        else
            resolution = (Int(ceil(resolution[1] * aspect / rescale)), resolution[2])
        end
    end

    # we need an integer resolution
    resolution = @. Int(round(resolution))

    # query vector field raster points
    rc, rv = vectorsample(grid, func; gridscale = ctx[:gridscale], rasterpoints = ((resolution[1] - 1) / 2, resolution[2] - 1), offset = ctx[:offset], xlimits = ex, ylimits = ey)
    qc, qv = quiverdata(rc, rv; vscale = ctx[:vscale], vnormalize = ctx[:vnormalize], vconstant = ctx[:vconstant])

    # construct canvas
    CanvasType = UnicodePlots.BrailleCanvas # should this be a changeable parameter ?
    canvas = CanvasType(
        resolution[2], resolution[1],               # number of rows and columns (characters)
        origin_y = ey[1], origin_x = ex[1],             # position in virtual space
        height = (ey[2] - ey[1]), width = (ex[2] - ex[1]); blend = false
    )

    # plot arrows
    narrows = size(qv, 2)
    vscale = ctx[:vscale] # vscale steers arrow thickness
    maxnorm = maximum(sqrt.(sum(qv .^ 2, dims = 1)))
    colormap = colorschemes[ctx[:colormap]]
    for a in 1:narrows
        # calculate angle of arrow
        angle = atan(qv[2, a], qv[1, a])
        anorm = sqrt(qv[1, a]^2 + qv[2, a]^2)
        scale = anorm / maxnorm
        uint_color = UnicodePlots.ansi_color(colormap[scale])
        char = select_arrow(angle, anorm / maxnorm, vscale)

        UnicodePlots.annotate!(canvas, qc[1, a], qc[2, a], char, uint_color, false)
    end

    # generate plot
    plt = UnicodePlots.Plot(
        canvas; title = ctx[:title], border = ctx[:border],
        xfact = (ex[2] - ex[1]) / (resolution[1] - 1),
        yfact = (ey[2] - ey[1]) / (resolution[2] - 1),
        xlabel = ctx[:xlabel],
        ylabel = ctx[:ylabel],
        xoffset = ex[1],
        yoffset = ey[1],
        compact_labels = false,
        labels = true
    )

    # add colormap
    plt.cmap.bar = ctx[:colorbar] == :none ? false : true
    plt.cmap.lim = (0, Float16(maxnorm))
    plt.cmap.callback = UnicodePlots.colormap_callback(ctx[:colormap])

    # corner coordinates
    ex = extrema(view(coords, 1, :))
    ey = extrema(view(coords, 2, :))
    UnicodePlots.label!(plt, :bl, UnicodePlots.nice_repr(ex[1], plt), UnicodePlots.ansi_color(UnicodePlots.BORDER_COLOR[]))
    UnicodePlots.label!(plt, :br, UnicodePlots.nice_repr(ex[2], plt), UnicodePlots.ansi_color(UnicodePlots.BORDER_COLOR[]))
    UnicodePlots.label!(plt, :l, 1, UnicodePlots.nice_repr(ey[2], plt), UnicodePlots.ansi_color(UnicodePlots.BORDER_COLOR[]))
    UnicodePlots.label!(plt, :l, resolution[2], UnicodePlots.nice_repr(ey[1], plt), UnicodePlots.ansi_color(UnicodePlots.BORDER_COLOR[]))

    ctx[:figure] = plt

    return reveal(ctx, TP)
end

streamplot!(ctx, TP::Type{UnicodePlotsType}, ::Type{Val{2}}, grid, func) = @warn "2D streamplot is not implemented for the UnicodePlots backend"


end # module
