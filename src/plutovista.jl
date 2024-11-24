using HypertextLiteral

function initialize!(p::GridVisualizer, ::Type{PlutoVistaType})
    PlutoVista = p.context[:Plotter]
    layout = p.context[:layout]
    figres = (p.context[:size][1] / layout[2], p.context[:size][2] / layout[1])
    for I in CartesianIndices(layout)
        ctx = p.subplots[I]
        ctx[:figure] = PlutoVista.PlutoVistaPlot(; resolution = figres)
        PlutoVista.backend!(
            ctx[:figure];
            backend = p.context[:backend],
            datadim = p.context[:dim],
        )
    end
    return nothing
end

function reveal(p::GridVisualizer, ::Type{PlutoVistaType})
    layout = p.context[:layout]
    figwidth = 0.95 * p.context[:size][1] / layout[2]
    l = layout[1] * layout[2]
    subplots = []
    for i in 1:layout[1]
        for j in 1:layout[2]
            push!(subplots, p.subplots[i, j])
        end
    end
    if l == 1
        return subplots[1][:figure]
    elseif l == 2
        return @htl(
            """<div><div style=" display: inline-block;">
            <div style= "width: $(figwidth)px;  display: inline-block;">$(subplots[1][:figure])</div>
            <div style= "width: $(figwidth)px;  display: inline-block;">$(subplots[2][:figure])</div>
            </div>"""
        )
    elseif l == 3
        return @htl(
            """<div><div style=" display: inline-block;">
            <div style= "width: $(figwidth)px;  display: inline-block;">$(subplots[1][:figure])</div>
            <div style= "width: $(figwidth)px;  display: inline-block;">$(subplots[2][:figure])</div>
            <div style= "width: $(figwidth)px;  display: inline-block;">$(subplots[3][:figure])</div>
            </div>"""
        )
    elseif l == 4
        return @htl(
            """<div><div style=" display: inline-block;">
            <div style= "width: $(figwidth)px;  display: inline-block;">$(subplots[1][:figure])</div>
            <div style= "width: $(figwidth)px;  display: inline-block;">$(subplots[2][:figure])</div>
            <div style= "width: $(figwidth)px;  display: inline-block;">$(subplots[3][:figure])</div>
            <div style= "width: $(figwidth)px;  display: inline-block;">$(subplots[4][:figure])</div>
            </div>"""
        )
    elseif l == 5
        return @htl(
            """<div><div style=" display: inline-block;">
            <div style= "width: $(figwidth)px;  display: inline-block;">$(subplots[1][:figure])</div>
            <div style= "width: $(figwidth)px;  display: inline-block;">$(subplots[2][:figure])</div>
            <div style= "width: $(figwidth)px;  display: inline-block;">$(subplots[3][:figure])</div>
            <div style= "width: $(figwidth)px;  display: inline-block;">$(subplots[4][:figure])</div>
            <div style= "width: $(figwidth)px;  display: inline-block;">$(subplots[5][:figure])</div>
            </div>"""
        )
    elseif l == 6
        return @htl(
            """<div><div style=" display: inline-block;">
            <div style= "width: $(figwidth)px;  display: inline-block;">$(subplots[1][:figure])</div>
            <div style= "width: $(figwidth)px;  display: inline-block;">$(subplots[2][:figure])</div>
            <div style= "width: $(figwidth)px;  display: inline-block;">$(subplots[3][:figure])</div>
            <div style= "width: $(figwidth)px;  display: inline-block;">$(subplots[4][:figure])</div>
            <div style= "width: $(figwidth)px;  display: inline-block;">$(subplots[5][:figure])</div>
            <div style= "width: $(figwidth)px;  display: inline-block;">$(subplots[6][:figure])</div>
            </div>"""
        )
    elseif l == 7
        return @htl(
            """<div><div style=" display: inline-block;">
            <div style= "width: $(figwidth)px;  display: inline-block;">$(subplots[1][:figure])</div>
            <div style= "width: $(figwidth)px;  display: inline-block;">$(subplots[2][:figure])</div>
            <div style= "width: $(figwidth)px;  display: inline-block;">$(subplots[3][:figure])</div>
            <div style= "width: $(figwidth)px;  display: inline-block;">$(subplots[4][:figure])</div>
            <div style= "width: $(figwidth)px;  display: inline-block;">$(subplots[5][:figure])</div>
            <div style= "width: $(figwidth)px;  display: inline-block;">$(subplots[6][:figure])</div>
            <div style= "width: $(figwidth)px;  display: inline-block;">$(subplots[7][:figure])</div>
            </div>"""
        )
    elseif l == 8
        return @htl(
            """<div><div style=" display: inline-block;">
            <div style= "width: $(figwidth)px;  display: inline-block;">$(subplots[1][:figure])</div>
            <div style= "width: $(figwidth)px;  display: inline-block;">$(subplots[2][:figure])</div>
            <div style= "width: $(figwidth)px;  display: inline-block;">$(subplots[3][:figure])</div>
            <div style= "width: $(figwidth)px;  display: inline-block;">$(subplots[4][:figure])</div>
            <div style= "width: $(figwidth)px;  display: inline-block;">$(subplots[5][:figure])</div>
            <div style= "width: $(figwidth)px;  display: inline-block;">$(subplots[6][:figure])</div>
            <div style= "width: $(figwidth)px;  display: inline-block;">$(subplots[7][:figure])</div>
            <div style= "width: $(figwidth)px;  display: inline-block;">$(subplots[8][:figure])</div>
            </div>"""
        )
    end
end

function reveal(ctx::SubVisualizer, TP::Type{PlutoVistaType})
    return if ctx[:show] || ctx[:reveal]
        reveal(ctx[:GridVisualizer], TP)
    end
end

function gridplot!(ctx, TP::Type{PlutoVistaType}, ::Type{Val{1}}, grid)
    PlutoVista = ctx[:Plotter]

    PlutoVista.backend!(ctx[:figure]; backend = ctx[:backend], datadim = 1)

    coord = grid[Coordinates]
    cellregions = cellcolors(grid, ctx[:cellcoloring])
    ncellregions = num_cellcolors(grid, ctx[:cellcoloring])
    cellnodes = grid[CellNodes]
    coord = grid[Coordinates] * ctx[:gridscale]
    bfacenodes = grid[BFaceNodes]
    bfaceregions = grid[BFaceRegions]
    nbfaceregions = grid[NumBFaceRegions]

    crflag = ones(Bool, ncellregions)
    brflag = ones(Bool, nbfaceregions)

    xmin = minimum(coord)
    xmax = maximum(coord)
    h = (xmax - xmin) / 20.0

    #    ax.set_aspect(ctx[:aspect])
    #    ax.get_yaxis().set_ticks([])
    #    ax.set_ylim(-5*h,xmax-xmin)
    cmap = region_cmap(ncellregions)

    for icell in 1:num_cells(grid)
        ireg = cellregions[icell]
        label = crflag[ireg] ? "c$(ireg)" : ""
        crflag[ireg] = false

        x1 = coord[1, cellnodes[1, icell]]
        x2 = coord[1, cellnodes[2, icell]]

        PlutoVista.plot!(
            ctx[:figure],
            [x1, x2],
            [0, 0];
            clear = false,
            linewidth = 3.0,
            color = rgbtuple(cmap[cellregions[icell]]),
            label = label,
        )
        PlutoVista.plot!(
            ctx[:figure],
            [x1, x1],
            [-h, h];
            clear = false,
            linewidth = ctx[:linewidth],
            color = :black,
        )
        PlutoVista.plot!(
            ctx[:figure],
            [x2, x2],
            [-h, h];
            clear = false,
            linewidth = ctx[:linewidth],
            color = :black,
        )
    end

    cmap = bregion_cmap(nbfaceregions)
    for ibface in 1:num_bfaces(grid)
        ireg = bfaceregions[ibface]
        if ireg > 0
            label = brflag[ireg] ? "b$(ireg)" : ""
            brflag[ireg] = false
            x1 = coord[1, bfacenodes[1, ibface]]
            PlutoVista.plot!(
                ctx[:figure],
                [x1, x1],
                [-2 * h, 2 * h];
                clear = false,
                linewidth = 3.0,
                color = rgbtuple(cmap[ireg]),
                label = label,
                legend = leglocs[ctx[:legend]],
            )
        end
    end

    return reveal(ctx, TP)
end

function scalarplot!(
        ctx,
        TP::Type{PlutoVistaType},
        ::Type{Val{1}},
        grids,
        parentgrid,
        funcs
    )
    nfuncs = length(funcs)

    PlutoVista = ctx[:Plotter]
    PlutoVista.backend!(ctx[:figure]; backend = ctx[:backend], datadim = 1)

    for ifunc in 1:nfuncs
        func = funcs[ifunc]
        grid = grids[ifunc]
        coord = grid[Coordinates] * ctx[:gridscale]

        if ifunc == 1
            PlutoVista.plot!(
                ctx[:figure],
                coord[1, :],
                func;
                color = ctx[:color],
                markertype = ctx[:markershape],
                markercount = length(func) ÷ ctx[:markevery],
                linestyle = ctx[:linestyle],
                xlabel = ctx[:xlabel],
                ylabel = ctx[:ylabel],
                label = ctx[:label],
                linewidth = ctx[:linewidth],
                legend = ctx[:legend],
                xlimits = ctx[:xlimits],
                limits = ctx[:limits],
                clear = ctx[:clear],
                title = ctx[:title],
                xscale = ctx[:xscale],
                yscale = ctx[:yscale],
            )
        else
            PlutoVista.plot!(
                ctx[:figure],
                coord[1, :],
                func;
                color = ctx[:color],
                markertype = ctx[:markershape],
                markercount = length(func) ÷ ctx[:markevery],
                linestyle = ctx[:linestyle],
                xlabel = ctx[:xlabel],
                ylabel = ctx[:ylabel],
                linewidth = ctx[:linewidth],
                xlimits = ctx[:xlimits],
                limits = ctx[:limits],
                clear = false,
                xscale = ctx[:xscale],
                yscale = ctx[:yscale],
            )
        end
    end
    return reveal(ctx, TP)
end

function gridplot!(ctx, TP::Type{PlutoVistaType}, ::Type{Val{2}}, grid)
    nregions = num_cellcolors(grid, ctx[:cellcoloring])
    nbregions = num_bfaceregions(grid)
    cmap = region_cmap(nregions)
    bcmap = bregion_cmap(nbregions)
    PlutoVista = ctx[:Plotter]
    pts = grid[Coordinates]
    tris = grid[CellNodes]
    markers = cellcolors(grid, ctx[:cellcoloring])
    edges = grid[BFaceNodes]
    edgemarkers = grid[BFaceRegions]

    if size(edges, 2) == 0
        edges = nothing
    end
    if size(edgemarkers, 1) == 0
        edgemarkers = nothing
        edgecolormap = nothing
    end

    PlutoVista.backend!(ctx[:figure]; backend = ctx[:backend], datadim = 2)
    PlutoVista.trimesh!(
        ctx[:figure],
        pts,
        tris;
        zoom = ctx[:zoom],
        markers = markers,
        colormap = cmap,
        gridscale = ctx[:gridscale],
        edges = edges,
        edgemarkers = edgemarkers,
        edgecolormap = bcmap,
        show_colorbar = ctx[:show_colorbar],
        aspect = ctx[:aspect],
        xlabel = ctx[:xlabel],
        ylabel = ctx[:ylabel],
        zlabel = ctx[:zlabel],
    )
    return reveal(ctx, TP)
end

function scalarplot!(
        ctx,
        TP::Type{PlutoVistaType},
        ::Type{Val{2}},
        grids,
        parentgrid,
        funcs
    )
    ngrids = length(grids)
    coords = [grid[Coordinates] for grid in grids]
    npoints = [num_nodes(grid) for grid in grids]
    cellnodes = [grid[CellNodes] for grid in grids]
    ncells = [num_cells(grid) for grid in grids]
    offsets = zeros(Int, ngrids)
    for i in 2:ngrids
        offsets[i] = offsets[i - 1] + npoints[i - 1]
    end

    allcoords = hcat(coords...)
    levels, crange, colorbarticks = isolevels(ctx, funcs)
    allcellnodes = Matrix{Int}(undef, 3, sum(ncells))
    k = 1
    for j in 1:ngrids
        for i in 1:ncells[j]
            allcellnodes[1, k] = cellnodes[j][1, i] + offsets[j]
            allcellnodes[2, k] = cellnodes[j][2, i] + offsets[j]
            allcellnodes[3, k] = cellnodes[j][3, i] + offsets[j]
            k = k + 1
        end
    end

    PlutoVista = ctx[:Plotter]
    PlutoVista.backend!(ctx[:figure]; backend = ctx[:backend], datadim = 2)
    PlutoVista.tricontour!(
        ctx[:figure],
        allcoords,
        allcellnodes,
        vcat(funcs...);
        colormap = ctx[:colormap],
        levels,
        colorbarticks = colorbarticks,
        limits = ctx[:limits],
        backend = ctx[:backend],
        zoom = ctx[:zoom],
        gridscale = ctx[:gridscale],
        aspect = ctx[:aspect],
        xlabel = ctx[:xlabel],
        ylabel = ctx[:ylabel],
        zlabel = ctx[:zlabel]
    )
    return reveal(ctx, TP)
end

function vectorplot!(ctx, TP::Type{PlutoVistaType}, ::Type{Val{2}}, grid, func)
    PlutoVista = ctx[:Plotter]
    PlutoVista.backend!(ctx[:figure]; backend = ctx[:backend], datadim = 2)

    rc, rv = vectorsample(grid, func; rasterpoints = ctx[:rasterpoints], offset = ctx[:offset], gridscale = ctx[:gridscale])
    qc, qv = quiverdata(rc, rv; vscale = ctx[:vscale], vnormalize = ctx[:vnormalize])

    PlutoVista.quiver2d!(ctx[:figure], qc, qv)
    return reveal(ctx, TP)
end

function streamplot!(ctx, TP::Type{PlutoVistaType}, ::Type{Val{2}}, grid, func) end

function gridplot!(ctx, TP::Type{PlutoVistaType}, ::Type{Val{3}}, grid)
    nbregions = num_bfaceregions(grid)
    nregions = num_cellcolors(grid, ctx[:cellcoloring])
    cmap = region_cmap(nregions)
    bcmap = bregion_cmap(nbregions)

    PlutoVista = ctx[:Plotter]
    pts = grid[Coordinates]
    tris = grid[CellNodes]
    faces = grid[BFaceNodes]
    markers = cellcolors(grid, ctx[:cellcoloring])
    facemarkers = grid[BFaceRegions]

    if size(faces, 2) == 0
        faces = nothing
    end
    if size(facemarkers, 1) == 0
        facemarkers = nothing
        facecolormap = nothing
    end


    PlutoVista.backend!(ctx[:figure]; backend = ctx[:backend], datadim = 3)
    PlutoVista.tetmesh!(
        ctx[:figure],
        pts,
        tris;
        xplanes = ctx[:xplanes],
        yplanes = ctx[:yplanes],
        zplanes = ctx[:zplanes],
        zoom = ctx[:zoom],
        markers = markers,
        colormap = cmap,
        faces = faces,
        facemarkers = facemarkers,
        gridscale = ctx[:gridscale],
        facecolormap = bcmap,
        show_colorbar = ctx[:show_colorbar],
        outlinealpha = ctx[:outlinealpha],
        xlabel = ctx[:xlabel],
        ylabel = ctx[:ylabel],
        zlabel = ctx[:zlabel],
    )
    return reveal(ctx, TP)
end

function scalarplot!(
        ctx,
        TP::Type{PlutoVistaType},
        ::Type{Val{3}},
        grids,
        parentgrid,
        funcs
    )
    PlutoVista = ctx[:Plotter]
    nbregions = num_bfaceregions(parentgrid)
    bcmap = bregion_cmap(nbregions)

    faces = parentgrid[BFaceNodes]
    facemarkers = parentgrid[BFaceRegions]
    facecolormap = bcmap

    if size(faces, 2) == 0
        faces = nothing
    end
    if size(facemarkers, 1) == 0
        facemarkers = nothing
        facecolormap = nothing
    end

    PlutoVista.backend!(ctx[:figure]; backend = ctx[:backend], datadim = 3)
    PlutoVista.tetcontour!(
        ctx[:figure],
        [grid[Coordinates] for grid in grids],
        [grid[CellNodes] for grid in grids],
        funcs;
        parentpts = parentgrid[Coordinates],
        levels = ctx[:levels],
        colormap = ctx[:colormap],
        xplanes = ctx[:xplanes],
        yplanes = ctx[:yplanes],
        zplanes = ctx[:zplanes],
        limits = ctx[:limits],
        faces = faces,
        facemarkers = facemarkers,
        facecolormap = facecolormap,
        outlinealpha = ctx[:outlinealpha],
        levelalpha = ctx[:levelalpha],
        zoom = ctx[:zoom],
        gridscale = ctx[:gridscale],
        aspect = ctx[:aspect],
        xlabel = ctx[:xlabel],
        ylabel = ctx[:ylabel],
        zlabel = ctx[:zlabel],
    )
    return reveal(ctx, TP)
end
