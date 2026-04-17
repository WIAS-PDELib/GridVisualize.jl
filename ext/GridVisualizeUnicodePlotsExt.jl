"""
    module GridVisualizeUnicodePlotsExt
    
Extension module for UnicodePlots.jl
"""
module GridVisualizeUnicodePlotsExt

import GridVisualize: initialize!, gridplot!, scalarplot!, bregion_cmap, region_cmap, reveal
using GridVisualize: UnicodePlotsType, GridVisualizer, SubVisualizer
using UnicodePlots: UnicodePlots
using ExtendableGrids: Coordinates, simplexgrid, ON_CELLS, ON_FACES, ON_EDGES, CellNodes, FaceNodes, BFaceNodes, CellGeometries, CellRegions, BFaceRegions, num_cells, num_nodes, local_celledgenodes, num_bfaceregions, num_cellregions, num_targets, interpolate!
using Colors: Colors, RGB, RGBA

initialize!(p, ::Type{UnicodePlotsType}) = nothing


function reveal(p::GridVisualizer, ::Type{UnicodePlotsType})
    layout = p.context[:layout]
    subplots = @views p.subplots[:]

    if layout == (1, 1)
        display(subplots[1][:figure])
    else
        if !isdefined(Main, :Term)
            @warn "A GridVisualizer with multiple UnicodePlots requires 'Term.jl' to be loaded: add Term.jl to your environment."
        else
            grid_plot = UnicodePlots.gridplot(map(subplot -> subplot[:figure], subplots), layout = p.context[:layout])
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


function region_legend!(canvas, title, x, y, colors)
    # legend by annotate!
    for (i, char) in enumerate(title)
        UnicodePlots.char_point!(canvas, x + i - 1, y, char, UInt32(0), false)
    end
    startx = x + length(title)
    for r in 1:length(colors)
        red, green, blue = UInt32.(colors[r])
        uint_color = (red << 16) | (green << 8) | blue
        reg_string = "$r "
        for char in reg_string
            startx += 1
            UnicodePlots.char_point!(canvas, startx, y, char, uint_color, false)
        end
    end
    return
end

function gridplot!(ctx, TP::Type{UnicodePlotsType}, ::Type{Val{2}}, grid)
    UnicodePlots = ctx[:Plotter]

    # find bounding box
    coords = grid[Coordinates]
    ex = extrema(view(coords, 1, :))
    ey = extrema(view(coords, 2, :))

    # line color for interior edges
    if typeof(ctx[:color]) <: RGB
        color = (
            Int(round(ctx[:color].r * 255)),
            Int(round(ctx[:color].g * 255)),
            Int(round(ctx[:color].b * 255)),
        )
    else
        color = ctx[:color]
    end

    # determine resolution (divided by 10, to reduce pixel count in the terminal)
    layout = ctx[:layout]
    resolution = ctx[:size] ./ 10 ./ (layout[2], layout[1])
    legend_space = 4
    aspect = ctx[:aspect] * resolution[1] / (resolution[1] + legend_space)

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

    # create UnicodePlots.Canvas
    padding = 0.1 * max(ex[2] - ex[1], ey[2] - ey[1])
    ex = (ex[1] - 2 * padding, ex[2] + 0.5 * padding)
    ey = (ey[1] - padding, ey[2] + padding)
    CanvasType = UnicodePlots.BrailleCanvas # should this be a changeable parameter ?
    canvas = CanvasType(
        resolution[2], resolution[1] + legend_space,            # number of rows and columns (characters)
        origin_y = ey[1], origin_x = ex[1],      # position in virtual space
        height = (ey[2] - ey[1]) / (resolution[1] / (resolution[1] + legend_space)), width = ex[2] - ex[1]; blend = false
    )

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
                color = color
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

    # color cell midpoints with cell regions color
    cellregions = grid[CellRegions]
    ncellregions = num_cellregions(grid)
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
    midpoint = [0.0, 0.0]
    for j in 1:ncells
        fill!(midpoint, 0.0)
        nvertices = num_targets(cellnodes, j)
        for k in 1:nvertices
            midpoint .+= coords[:, cellnodes[k, j]]
        end
        midpoint ./= nvertices
        r = cellregions[j]
        UnicodePlots.points!(
            canvas,
            midpoint[1], midpoint[2];
            color = colors[r]
        )
    end

    # plot boundary faces with bregion_cmap colors
    nbregions = num_bfaceregions(grid)
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

    region_legend!(canvas, "cell regions: ", 2, 1, colors)
    region_legend!(canvas, "bface regions:", 2, 2, bcolors)

    # corner coordinates
    ex = extrema(view(coords, 1, :))
    ey = extrema(view(coords, 2, :))
    UnicodePlots.annotate!(canvas, ex[1], ey[1], "$(ex[1])", UInt32(0), false; valign = :top)
    UnicodePlots.annotate!(canvas, ex[2], ey[1], "$(ex[2])", UInt32(0), false; valign = :top, halign = :right)
    UnicodePlots.annotate!(canvas, ex[1] - 1.5 * padding, ey[1], "$(ey[1])", UInt32(0), false; halign = :left)
    UnicodePlots.annotate!(canvas, ex[1] - 1.5 * padding, ey[2], "$(ey[2])", UInt32(0), false; halign = :left)

    # plot
    ctx[:figure] = UnicodePlots.Plot(canvas; title = ctx[:title])
    return reveal(ctx, TP)
end


function gridplot!(ctx, TP::Type{UnicodePlotsType}, ::Type{Val{1}}, grid)
    UnicodePlots = ctx[:Plotter]

    # find bounding box
    coords = grid[Coordinates]
    ex = extrema(view(coords, 1, :))

    # line color for interior edges
    if typeof(ctx[:color]) <: RGB
        color = (
            Int(round(ctx[:color].r * 255)),
            Int(round(ctx[:color].g * 255)),
            Int(round(ctx[:color].b * 255)),
        )
    else
        color = ctx[:color]
    end

    # determine resolution (divided by 10, to reduce pixel count in the terminal)
    layout = ctx[:layout]
    resolution = (Int(round(ctx[:size][1] / 10 / layout[2])), 5)

    # create UnicodePlots.Canvas
    legend_space = 5
    padding = 0.05 * (ex[2] - ex[1])
    ex = (ex[1] - padding, ex[2] + padding)
    CanvasType = UnicodePlots.BrailleCanvas # should this be a changeable parameter ?
    canvas = CanvasType(
        resolution[2], resolution[1] + legend_space,               # number of rows and columns (characters)
        origin_y = 0, origin_x = ex[1],             # position in virtual space
        height = 1, width = ex[2] - ex[1]; blend = false
    )

    # plot all edges in the grid
    cellregions = grid[CellRegions]
    ncellregions = num_cellregions(grid)
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
    for j in 1:ncells
        cen = local_celledgenodes(cellgeoms[j])
        r = cellregions[j]
        for k in 1:size(cen, 2)
            UnicodePlots.lines!(
                canvas,
                coords[1, cellnodes[cen[1, k], j]], 0.3,
                coords[1, cellnodes[cen[2, k], j]], 0.3;
                color = colors[r]
            )
        end
    end
    for j in 1:nnodes
        UnicodePlots.annotate!(canvas, coords[1, j], 0.4, "•", UInt32(0), false)
    end

    # plot boundary nodes with bregion_cmap colors
    nbregions = num_bfaceregions(grid)
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
        UnicodePlots.annotate!(canvas, coords[1, bfacenodes[1, j]], 0.4, "•", UInt32(uint_color), false)
    end

    region_legend!(canvas, "cell regions: ", 2, 1, colors)
    region_legend!(canvas, "bface regions:", 2, 2, bcolors)


    ex = extrema(view(coords, 1, :))
    UnicodePlots.annotate!(canvas, 0, 0.1, "$(ex[1])", UInt32(0), false)
    UnicodePlots.annotate!(canvas, ex[2], 0.1, "$(ex[2])", UInt32(0), false)

    # plot
    ctx[:figure] = UnicodePlots.Plot(canvas; title = ctx[:title])

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
    resolution = @. Int(round(ctx[:size] ./ 10 ./ (layout[2], layout[1]))) # reduce pixel count in the terminal
    ylim = ctx[:limits]

    if ylim[1] > ylim[2]
        # try to find limits automatically
        ylim = (minimum([minimum(func) for func in funcs]), maximum([maximum(func) for func in funcs]))
    end

    plt = ctx[:clear] ? nothing : ctx[:figure]
    for ifunc in 1:nfuncs
        func = funcs[ifunc]
        grid = grids[ifunc]
        coord = grid[Coordinates] * ctx[:gridscale]
        if ifunc == 1
            plt = UnicodePlots.lineplot(coord[1, :], func; ylim, xlabel = "x", name = ctx[:label], height = resolution[2], width = resolution[1])
        else
            UnicodePlots.lineplot!(plt, coord[1, :], func; name = ctx[:label])
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
    resolution = ctx[:size] ./ 10 ./ (layout[2], layout[1]) # reduce pixel count in the terminal
    ylim = ctx[:limits]
    colormap = ctx[:colormap]

    if ylim[1] > ylim[2]
        # try to find limits automatically
        ylim = (minimum([minimum(func) for func in funcs]), maximum([maximum(func) for func in funcs]))
    end

    coords = grids[1][Coordinates]
    ex = extrema(view(coords, 1, :))
    ey = extrema(view(coords, 2, :))

    if (true) # auto scale feature, do we want this?
        wx = ex[2] - ex[1]
        wy = ey[2] - ey[1]
        rescale = wx / wy * (resolution[1] / (resolution[2]))
        if rescale > 1
            resolution = (resolution[1], resolution[2] / rescale)
        else
            resolution = (resolution[1] / rescale, resolution[2])
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

    ctx[:figure] = UnicodePlots.heatmap(
        reshape(I, (resolution[1], resolution[2]))',
        xlabel = "x",
        ylabel = "y",
        xfact = (ex[2] - ex[1]) / (resolution[1] - 1),
        yfact = (ey[2] - ey[1]) / (resolution[2] - 1),
        xoffset = ex[1],
        yoffset = ey[1],
        title = ctx[:title],
        colormap = colormap,
    )

    return reveal(ctx, TP)
end

end # module
