function initialize!(p::GridVisualizer, ::Type{VTKViewType})
    pctx = p.context
    VTKView = pctx[:Plotter]
    frame = VTKView.StaticFrame()
    pctx[:frame] = frame
    if pctx[:clear]
        VTKView.clear!(frame)
    end
    layout = pctx[:layout]
    tlayout = (layout[2], layout[1])
    for I in CartesianIndices(layout)
        ctx = p.subplots[I]
        ctx[:frame] = frame
    end
    VTKView.layout!(frame, tlayout...)
    VTKView.size!(frame, pctx[:size]...)
    return pctx
end

function save(fname, p, ::Type{VTKViewType})
    VTKView = p.context[:Plotter]
    base, ext = splitext(fname)
    if ext != ".png"
        error("VTKView can only save png files")
    end
    return VTKView.writepng(p.context[:frame], fname)
end

function save(fname, scene, VTKView, ::Type{VTKViewType})
    base, ext = splitext(fname)
    if ext != ".png"
        error("VTKView can only save png files")
    end
    return VTKView.writepng(scene, fname)
end

function reveal(p::GridVisualizer, ::Type{VTKViewType})
    VTKView = p.Plotter
    VTKView.display(p.context[:frame])
    return p.context[:frame]
end

function reveal(ctx::SubVisualizer, TP::Type{VTKViewType})
    return if ctx[:show] || ctx[:reveal]
        reveal(ctx[:GridVisualizer], TP)
    end
end

function gridplot!(ctx, TP::Type{VTKViewType}, grid)
    VTKView = ctx[:Plotter]
    frame = ctx[:frame]
    if !haskey(ctx, :dataset)
        dataset = VTKView.DataSet()
        ctx[:dataset] = dataset
    end
    VTKView.simplexgrid!(ctx[:dataset], grid[Coordinates], grid[CellNodes])
    VTKView.boundarygrid!(ctx[:dataset], grid[BFaceNodes])
    VTKView.boundarymarker!(ctx[:dataset], grid[BFaceRegions])
    VTKView.cellmarker!(ctx[:dataset], cellcolors(grid, ctx[:cellcoloring]))
    if !haskey(ctx, :gridview)
        ctx[:gridview] = VTKView.GridView()
        VTKView.data!(ctx[:gridview], ctx[:dataset])
        VTKView.addview!(frame, ctx[:gridview], ctx[:iplot]...)
    end
    return reveal(ctx, TP)
end

gridplot!(ctx, T::Type{VTKViewType}, ::Type{Val{2}}, grid) = gridplot!(ctx, T, grid)
gridplot!(ctx, T::Type{VTKViewType}, ::Type{Val{3}}, grid) = gridplot!(ctx, T, grid)

function scalarplot!(ctx, TP::Type{VTKViewType}, grids, parentgrid, funcs)
    grid = parentgrid
    func = funcs[1]
    VTKView = ctx[:Plotter]
    frame = ctx[:frame]
    if !haskey(ctx, :dataset)
        ctx[:dataset] = VTKView.DataSet()
    end
    if !haskey(ctx, :grid)
        ctx[:grid] = grid
        VTKView.simplexgrid!(ctx[:dataset], grid[Coordinates], grid[CellNodes])
    end
    if !haskey(ctx, :scalarview)
        ctx[:scalarview] = VTKView.ScalarView()
        VTKView.addview!(frame, ctx[:scalarview], ctx[:iplot])
        VTKView.data!(ctx[:scalarview], ctx[:dataset], ctx[:label])
    end
    if !seemingly_equal(grid, ctx[:grid])
        VTKView.simplexgrid!(ctx[:dataset], grid[Coordinates], grid[CellNodes])
    end
    VTKView.pointscalar!(ctx[:dataset], func, ctx[:label])
    if dim_space(grid) == 3
        VTKView.isolevels!(ctx[:scalarview], [ctx[:flevel]])
        VTKView.show_isosurfaces!(ctx[:scalarview], true)
    end
    return reveal(ctx, TP)
end

function scalarplot!(ctx, T::Type{VTKViewType}, ::Type{Val{2}}, grids, parentgrid, funcs)
    return scalarplot!(ctx, T, grids, parentgrid, funcs)
end
function scalarplot!(ctx, T::Type{VTKViewType}, ::Type{Val{3}}, grids, parentgrid, funcs)
    return scalarplot!(ctx, T, grids, parentgrid, funcs)
end
gridplot!(ctx, T::Type{VTKViewType}, ::Type{Val{1}}, grid) = nothing

function scalarplot!(ctx, TP::Type{VTKViewType}, ::Type{Val{1}}, grids, parentgrid, funcs)
    grid = parentgrid
    func = funcs[1]
    VTKView = ctx[:Plotter]
    frame = ctx[:frame]
    if !haskey(ctx, :plot)
        ctx[:plot] = VTKView.XYPlot()
        VTKView.addview!(frame, ctx[:plot], ctx[:iplot])
        VTKView.xrange!(ctx[:plot], ctx[:xlimits]...)
        VTKView.yrange!(ctx[:plot], ctx[:limits]...)
        VTKView.linewidth!(ctx[:plot], 1)
    end
    if ctx[:clear]
        VTKView.clear!(ctx[:plot])
    end
    plot = ctx[:plot]
    VTKView.plotcolor!(plot, rgbtuple(ctx[:color])...)

    VTKView.title!(plot, ctx[:title])

    if ctx[:label] != ""
        VTKView.plotlegend!(plot, ctx[:label])
    end

    VTKView.addplot!(plot, collect(grid[Coordinates][1, :]), collect(func))
    return reveal(ctx, TP)
end

vectorplot!(ctx, TP::Type{VTKViewType}, ::Type{Val{1}}, grid, func) = nothing
vectorplot!(ctx, TP::Type{VTKViewType}, ::Type{Val{2}}, grid, func) = nothing
vectorplot!(ctx, TP::Type{VTKViewType}, ::Type{Val{3}}, grid, func) = nothing

streamplot!(ctx, TP::Type{VTKViewType}, ::Type{Val{1}}, grid, func) = nothing
streamplot!(ctx, TP::Type{VTKViewType}, ::Type{Val{2}}, grid, func) = nothing
streamplot!(ctx, TP::Type{VTKViewType}, ::Type{Val{3}}, grid, func) = nothing
