function initialize!(p, ::Type{PyPlotType})
    PyPlot = p.context[:Plotter]
    PyPlot.PyObject(PyPlot.axes3D) # see https://github.com/JuliaPy/PyPlot.jl/issues/351
    PyPlot.rc("font"; size = p.context[:fontsize])
    if !haskey(p.context, :figure)
        res = p.context[:size]
        if !isdefined(Main, :PlutoRunner)
            p.context[:figure] = PyPlot.figure(p.context[:fignumber]; dpi = 50, figsize = res ./ 50)
        else
            p.context[:figure] = PyPlot.figure(p.context[:fignumber]; dpi = 100, figsize = res ./ 100)
        end
        #p.context[:figure].set_size_inches(res[1] / 100, res[2] / 100, forward = true)
        for ctx in p.subplots
            ctx[:figure] = p.context[:figure]
        end
    end
    return if p.context[:clear]
        p.context[:figure].clf()
        p.context[:revealed] = false
    end
end

function save(fname, p, ::Type{PyPlotType})
    return p.context[:figure].savefig(fname)
end

function save(fname, scene, PyPlot, ::Type{PyPlotType})
    return isnothing(scene) ? nothing : scene.savefig(fname)
end

function reveal(p::GridVisualizer, ::Type{PyPlotType})
    p.context[:revealed] = true
    p.Plotter.tight_layout()
    if !(isdefined(Main, :PlutoRunner)) && isinteractive()
        p.Plotter.pause(1.0e-10)
        p.Plotter.draw()
        sleep(1.0e-3)
    end
    return p.context[:figure]
end

function reveal(ctx::SubVisualizer, TP::Type{PyPlotType})
    yield()
    if ctx[:show] || ctx[:reveal]
        reveal(ctx[:GridVisualizer], TP)
    end
    return ctx[:GridVisualizer].Plotter.tight_layout()
end

#translate Julia attribute symbols to pyplot-speak
const mshapes = Dict(
    :dtriangle => "v",
    :utriangle => "^",
    :rtriangle => ">",
    :ltriangle => "^",
    :circle => "o",
    :square => "s",
    :cross => "+",
    :+ => "+",
    :xcross => "x",
    :x => "x",
    :diamond => "D",
    :star5 => "*",
    :pentagon => "p",
    :hexagon => "h"
)

const lstyles = Dict(
    :solid => "-",
    :dot => "dotted",
    :dash => "--",
    :dashdot => "-.",
    :dashdotdot => (0, (3, 1, 1, 1))
)

const leglocs = Dict(
    :none => "",
    :best => "best",
    :lt => "upper left",
    :ct => "upper center",
    :rt => "upper right",
    :lc => "center left",
    :cc => "center center",
    :rc => "center right",
    :lb => "lower left",
    :cb => "lower center",
    :rb => "lower right"
)

"""
$(SIGNATURES)
Return tridata to be splatted to PyPlot calls
"""
function tridata(grid::ExtendableGrid, gridscale)
    coord = grid[Coordinates] * gridscale
    cellnodes = Matrix(grid[CellNodes])
    return coord[1, :], coord[2, :], transpose(cellnodes .- 1)
end

function tridata(grids, gridscale)
    ngrids = length(grids)
    coords = [grid[Coordinates] * gridscale for grid in grids]
    npoints = [num_nodes(grid) for grid in grids]
    cellnodes = [grid[CellNodes] for grid in grids]
    ncells = [num_cells(grid) for grid in grids]
    offsets = zeros(Int, ngrids)
    for i in 2:ngrids
        offsets[i] = offsets[i - 1] + npoints[i - 1]
    end

    allcoords = hcat(coords...)

    # transpose and subtract 1 !
    allcellnodes = Matrix{Int}(undef, sum(ncells), 3)
    k = 1
    for j in 1:ngrids
        for i in 1:ncells[j]
            allcellnodes[k, 1] = cellnodes[j][1, i] + offsets[j] - 1
            allcellnodes[k, 2] = cellnodes[j][2, i] + offsets[j] - 1
            allcellnodes[k, 3] = cellnodes[j][3, i] + offsets[j] - 1
            k = k + 1
        end
    end
    return allcoords[1, :], allcoords[2, :], allcellnodes
end

# Interfaces to Colors/Colorschemes
plaincolormap(ctx) = colorschemes[ctx[:colormap]].colors

### 1D grid
function gridplot!(ctx, TP::Type{PyPlotType}, ::Type{Val{1}}, grid)
    PyPlot = ctx[:Plotter]

    if !haskey(ctx, :ax)
        ctx[:ax] = ctx[:figure].add_subplot(ctx[:layout]..., ctx[:iplot])
    end
    if ctx[:clear]
        ctx[:ax].cla()
    end
    ax = ctx[:ax]
    fig = ctx[:figure]

    cellregions = cellcolors(grid, ctx[:cellcoloring])
    ncellregions = num_cellcolors(grid, ctx[:cellcoloring])
    cellnodes = grid[CellNodes]
    coord = grid[Coordinates]
    bfacenodes = grid[BFaceNodes]
    bfaceregions = grid[BFaceRegions]
    nbfaceregions = grid[NumBFaceRegions]

    crflag = ones(Bool, ncellregions)
    brflag = ones(Bool, nbfaceregions)

    xmin = minimum(coord)
    xmax = maximum(coord)
    h = (xmax - xmin) / 20.0
    #    ax.set_aspect(ctx[:aspect])
    ax.grid(true)
    ax.get_yaxis().set_ticks([])
    ax.set_ylim(-5 * h, xmax - xmin)
    cmap = region_cmap(ncellregions)
    gridscale = ctx[:gridscale]

    for icell in 1:num_cells(grid)
        ireg = cellregions[icell]
        label = crflag[ireg] ? "c$(ireg)" : ""
        crflag[ireg] = false

        x1 = coord[1, cellnodes[1, icell]] * gridscale
        x2 = coord[1, cellnodes[2, icell]] * gridscale
        ax.plot(
            [x1, x2],
            [0, 0];
            linewidth = 3.0,
            color = rgbtuple(cmap[cellregions[icell]]),
            label = label,
        )
        ax.plot([x1, x1], [-h, h]; linewidth = ctx[:linewidth], color = "k", label = "")
        ax.plot([x2, x2], [-h, h]; linewidth = ctx[:linewidth], color = "k", label = "")
    end

    cmap = bregion_cmap(max(nbfaceregions,2))
    for ibface in 1:num_bfaces(grid)
        ireg = bfaceregions[ibface]
        if ireg > 0
            label = brflag[ireg] ? "b$(ireg)" : ""
            brflag[ireg] = false
            x1 = coord[1, bfacenodes[1, ibface]] * ctx[:gridscale]
            ax.plot(
                [x1, x1],
                [-2 * h, 2 * h];
                linewidth = 3.0,
                color = rgbtuple(cmap[ireg]),
                label = label,
            )
        end
    end
    if ctx[:legend] != :none
        ax.legend(; loc = leglocs[ctx[:legend]], ncol = 5)
    end
    return reveal(ctx, TP)
end

### 2D grid
function gridplot!(ctx, TP::Type{PyPlotType}, ::Type{Val{2}}, grid)
    PyPlot = ctx[:Plotter]
    if !haskey(ctx, :ax)
        ctx[:ax] = ctx[:figure].add_subplot(ctx[:layout]..., ctx[:iplot])
    end
    if ctx[:clear]
        ctx[:ax].cla()
        ctx[:ax].set_title(ctx[:title])

        xlimits = ctx[:xlimits]
        ylimits = ctx[:ylimits]
        if xlimits[1] < xlimits[2]
            ctx[:ax].set_xlim(xlimits...)
        end
        if ylimits[1] < ylimits[2]
            ctx[:ax].set_ylim(ylimits...)
        end
    end
    ax = ctx[:ax]
    fig = ctx[:figure]
    cellregions = cellcolors(grid, ctx[:cellcoloring])
    ncellregions = num_cellcolors(grid, ctx[:cellcoloring])
    cellnodes = grid[CellNodes]
    nbfaceregions = grid[NumBFaceRegions]
    if nbfaceregions > 0
        bfacenodes = grid[BFaceNodes]
        bfaceregions = grid[BFaceRegions]
    end

    crflag = ones(Bool, ncellregions)
    brflag = ones(Bool, nbfaceregions)
    ax.set_aspect(ctx[:aspect])
    tridat = tridata(grid, ctx[:gridscale])
    cmap = region_cmap(max(ncellregions, 2))
    cdata = ax.tripcolor(
        tridat...;
        facecolors = cellcolors(grid, ctx[:cellcoloring]),
        cmap = PyPlot.ColorMap(cmap, length(cmap)),
        vmin = 1.0,
        vmax = length(cmap),
    )
    if ctx[:show_colorbar]
        if ctx[:colorbar] == :horizontal
            cbar = fig.colorbar(
                cdata;
                ax = ax,
                ticks = collect(1:length(cmap)),
                orientation = "horizontal",
            )
        end

        if ctx[:colorbar] == :vertical
            cbar = fig.colorbar(
                cdata;
                ax = ax,
                ticks = collect(1:length(cmap)),
                orientation = "vertical",
            )
        end
    end

    ax.triplot(tridat...; color = "k", linewidth = ctx[:linewidth])

    if nbfaceregions > 0
        gridscale = ctx[:gridscale]
        coord = grid[Coordinates]
        cmap = bregion_cmap(max(nbfaceregions, 2))
        # see https://gist.github.com/gizmaa/7214002
        c1 = [coord[:, bfacenodes[1, i]] for i in 1:num_sources(bfacenodes)] * gridscale
        c2 = [coord[:, bfacenodes[2, i]] for i in 1:num_sources(bfacenodes)] * gridscale
        rgb = [rgbtuple(cmap[bfaceregions[i]]) for i in 1:length(bfaceregions)]
        ax.add_collection(
            PyPlot.matplotlib.collections.LineCollection(
                collect(zip(c1, c2));
                colors = rgb,
                linewidth = 3,
            )
        )
        for i in 1:nbfaceregions
            ax.plot(coord[1, 1:1] * gridscale, coord[2, 1:1] * gridscale; label = "$(i)", color = rgbtuple(cmap[i]))
        end
    end
    if ctx[:legend] != :none
        ax.legend(; loc = leglocs[ctx[:legend]])
    end
    return reveal(ctx, TP)
end

### 3D Grid
function gridplot!(ctx, TP::Type{PyPlotType}, ::Type{Val{3}}, grid)
    # See https://jakevdp.github.io/PythonDataScienceHandbook/04.12-three-dimensional-plotting.html

    PyPlot = ctx[:Plotter]
    if !haskey(ctx, :ax)
        ctx[:ax] = ctx[:figure].add_subplot(ctx[:layout]..., ctx[:iplot]; projection = "3d")
    end

    ax = ctx[:ax]

    if ctx[:clear]
        ctx[:ax].cla()
    end
    fig = ctx[:figure]

    nregions = num_cellregions(grid)
    nbregions = num_bfaceregions(grid)
    gridscale = ctx[:gridscale]
    xyzmin = zeros(3)
    xyzmax = ones(3)
    coord = grid[Coordinates] * gridscale
    @views for idim in 1:3
        xyzmin[idim] = minimum(coord[idim, :])
        xyzmax[idim] = maximum(coord[idim, :])
    end

    ax.set_xlim3d(xyzmin[1], xyzmax[1])
    ax.set_ylim3d(xyzmin[2], xyzmax[2])
    ax.set_zlim3d(xyzmin[3], xyzmax[3])
    ax.view_init(ctx[:elev], ctx[:azim])

    cmap = region_cmap(max(nregions, 2))
    bcmap = bregion_cmap(max(nbregions, 2))

    xyzcut = [ctx[:xplanes][1], ctx[:yplanes][1], ctx[:zplanes][1]]

    if ctx[:interior]
        regpoints0, regfacets0 = extract_visible_cells3D(
            grid, xyzcut; gridscale = ctx[:gridscale],
            cellcoloring = ctx[:cellcoloring],
            primepoints = hcat(xyzmin, xyzmax)
        )
        regfacets = [
            reshape(reinterpret(Int32, regfacets0[i]), (3, length(regfacets0[i]))) for
                i in 1:nregions
        ]
        regpoints = [
            reshape(reinterpret(Float32, regpoints0[i]), (3, length(regpoints0[i]))) for
                i in 1:nregions
        ]

        for ireg in 1:nregions
            if size(regfacets[ireg], 2) > 0
                ax.plot_trisurf(
                    regpoints[ireg][1, :],
                    regpoints[ireg][2, :],
                    transpose(regfacets[ireg] .- 1),
                    regpoints[ireg][3, :];
                    color = rgbtuple(cmap[ireg]),
                    edgecolors = :black,
                    linewidth = 0.5,
                )
            end
        end
    end

    bregpoints0, bregfacets0 = extract_visible_bfaces3D(
        grid, xyzcut; gridscale = ctx[:gridscale],
        primepoints = hcat(xyzmin, xyzmax)
    )
    bregfacets = [
        reshape(reinterpret(Int32, bregfacets0[i]), (3, length(bregfacets0[i]))) for
            i in 1:nbregions
    ]
    bregpoints = [
        reshape(reinterpret(Float32, bregpoints0[i]), (3, length(bregpoints0[i]))) for
            i in 1:nbregions
    ]
    for ireg in 1:nbregions
        if size(bregfacets[ireg], 2) > 0
            ax.plot_trisurf(
                bregpoints[ireg][1, :],
                bregpoints[ireg][2, :],
                transpose(bregfacets[ireg] .- 1),
                bregpoints[ireg][3, :];
                color = rgbtuple(bcmap[ireg]),
                edgecolors = :black,
                linewidth = 0.5,
            )
        end
    end

    if ctx[:legend] != :none
        ax.legend(; loc = leglocs[ctx[:legend]])
    end
    return reveal(ctx, TP)
end

### 1D Function
function scalarplot!(ctx, TP::Type{PyPlotType}, ::Type{Val{1}}, grids, parentgrid, funcs)
    PyPlot = ctx[:Plotter]
    nfuncs = length(funcs)

    function set_limits_grid_title()
        ax = ctx[:ax]
        xlimits = ctx[:xlimits]
        ylimits = ctx[:limits]

        if xlimits[1] < xlimits[2]
            ax.set_xlim(xlimits...)
        end
        if ylimits[1] < ylimits[2]
            ax.set_ylim(ylimits...)
        end
        ax.set_title(ctx[:title])
        return ax.grid()
    end

    if !haskey(ctx, :ax)
        ctx[:ax] = ctx[:figure].add_subplot(ctx[:layout]..., ctx[:iplot])
        set_limits_grid_title()
    end
    if ctx[:clear]
        ctx[:ax].cla()
        set_limits_grid_title()
    end
    ax = ctx[:ax]
    fig = ctx[:figure]

    pplot = ax.plot
    if ctx[:xscale] == :log
        if ctx[:yscale] == :log
            pplot = ax.loglog
        else
            pplot = ax.semilogx
        end
    end
    if ctx[:yscale] == :log
        if ctx[:xscale] == :log
            pplot = ax.loglog
        else
            pplot = ax.semilogy
        end
    end
    gridscale = ctx[:gridscale]
    if ctx[:cellwise] # not checked,  outdated
        for icell in 1:num_cells(grid)
            i1 = cellnodes[1, icell]
            i2 = cellnodes[2, icell]
            x1 = coord[1, i1]
            x2 = coord[1, i2]
            if icell == 1
                ax.plot(
                    [x1, x2],
                    [func[i1], func[i2]];
                    color = rgbtuple(ctx[:color]),
                    label = ctx[:label],
                )
            else
                ax.plot([x1, x2], [func[i1], func[i2]]; color = rgbtuple(ctx[:color]))
            end
        end
    else
        if ctx[:markershape] == :none
            for ifunc in 1:nfuncs
                func = funcs[ifunc]
                coord = grids[ifunc][Coordinates] * gridscale
                if ctx[:label] !== "" && ifunc == 1
                    pplot(
                        coord[1, :],
                        func;
                        linestyle = lstyles[ctx[:linestyle]],
                        color = rgbtuple(ctx[:color]),
                        linewidth = ctx[:linewidth],
                        label = ctx[:label],
                    )
                else
                    pplot(
                        coord[1, :],
                        func;
                        linestyle = lstyles[ctx[:linestyle]],
                        linewidth = ctx[:linewidth],
                        color = rgbtuple(ctx[:color]),
                    )
                end
            end
        else
            for ifunc in 1:nfuncs
                func = funcs[ifunc]
                coord = grids[ifunc][Coordinates] * gridscale
                if ctx[:label] !== "" && ifunc == 1
                    pplot(
                        coord[1, :],
                        func;
                        linestyle = lstyles[ctx[:linestyle]],
                        color = rgbtuple(ctx[:color]),
                        label = ctx[:label],
                        marker = mshapes[ctx[:markershape]],
                        markevery = ctx[:markevery],
                        markersize = ctx[:markersize],
                        linewidth = ctx[:linewidth],
                    )
                else
                    pplot(
                        coord[1, :],
                        func;
                        linestyle = lstyles[ctx[:linestyle]],
                        color = rgbtuple(ctx[:color]),
                        marker = mshapes[ctx[:markershape]],
                        markevery = ctx[:markevery],
                        markersize = ctx[:markersize],
                        linewidth = ctx[:linewidth],
                    )
                end
            end
        end
        # points=[Point2f(coord[1,i],func[i]) for i=1:length(func)]
        # Hard to get this robust, as we need to get axislimits
        # mpoints=markerpoints(points,ctx[:markers],Diagonal([1,1]))
        # ampoints=reshape(reinterpret(Float32,mpoints),(2,length(mpoints)))
        # ax.scatter(ampoints[1,:], ampoints[2,:],color=ctx[:color],label="")

    end
    ax.grid(true)
    ax.set_xlabel(ctx[:xlabel])
    ax.set_ylabel(ctx[:ylabel])
    if ctx[:legend] != :none
        ax.legend(; loc = leglocs[ctx[:legend]])
    end

    return reveal(ctx, TP)
end

### 2D Function
function scalarplot!(ctx, TP::Type{PyPlotType}, ::Type{Val{2}}, grids, parentgrid, funcs)
    PyPlot = ctx[:Plotter]
    if !haskey(ctx, :ax)
        ctx[:ax] = ctx[:figure].add_subplot(ctx[:layout]..., ctx[:iplot])
    end
    if ctx[:clear]
        if haskey(ctx, :cbar)
            ctx[:cbar].remove()
        end
        ctx[:ax].remove()
        ctx[:ax] = ctx[:figure].add_subplot(ctx[:layout]..., ctx[:iplot])

        xlimits = ctx[:xlimits]
        ylimits = ctx[:ylimits]
        if xlimits[1] < xlimits[2]
            ctx[:ax].set_xlim(xlimits...)
        end
        if ylimits[1] < ylimits[2]
            ctx[:ax].set_ylim(ylimits...)
        end
    end

    ax = ctx[:ax]
    fig = ctx[:figure]
    ax.set_aspect(ctx[:aspect])
    ax.set_title(ctx[:title])

    levels, crange, colorbarticks = isolevels(ctx, funcs)
    eps = 1.0e-5
    if crange[1] == crange[2]
        eps = 1.0e-5
    else
        eps = (crange[2] - crange[1]) * 1.0e-15
    end

    colorlevels = range(crange[1] - eps, crange[2] + eps; length = ctx[:colorlevels])

    #    if !haskey(ctx, :grid) || !seemingly_equal(ctx[:grid], grid)
    #        ctx[:grid] = grids
    #        ctx[:tridata] = tridata(grids)
    #    end

    tdat = tridata(grids, ctx[:gridscale])
    func = vcat(funcs...)
    cnt = ax.tricontourf(
        tdat...,
        func;
        levels = colorlevels,
        cmap = PyPlot.ColorMap(plaincolormap(ctx)),
    )

    for c in cnt.collections
        c.set_edgecolor("face")
    end

    ax.tricontour(tdat..., func; colors = "k", levels = levels)

    if ctx[:colorbar] == :horizontal
        ctx[:cbar] = fig.colorbar(
            cnt;
            ax = ax,
            ticks = colorbarticks,
            boundaries = colorlevels,
            orientation = "horizontal",
        )
    end
    if ctx[:colorbar] == :vertical
        ctx[:cbar] = fig.colorbar(
            cnt;
            ax = ax,
            ticks = colorbarticks,
            boundaries = colorlevels,
            orientation = "vertical",
        )
    end

    ax.set_xlabel(ctx[:xlabel])

    ax.set_ylabel(ctx[:ylabel])

    return reveal(ctx, TP)
end

function scalarplot!(ctx, TP::Type{PyPlotType}, ::Type{Val{3}}, grids, parentgrid, funcs)
    PyPlot = ctx[:Plotter]
    if !haskey(ctx, :ax)
        ctx[:ax] = ctx[:figure].add_subplot(ctx[:layout]..., ctx[:iplot]; projection = "3d")
    end
    ax = ctx[:ax]
    fig = ctx[:figure]
    griscale = ctx[:gridscale]
    xyzmin = zeros(3)
    xyzmax = ones(3)
    coord = parentgrid[Coordinates]
    @views for idim in 1:3
        xyzmin[idim] = minimum(coord[idim, :]) * ctx[:gridscale]
        xyzmax[idim] = maximum(coord[idim, :]) * ctx[:gridscale]
    end
    xyzcut = [ctx[:xplanes], ctx[:yplanes], ctx[:zplanes]] * ctx[:gridscale]
    levels, crange, colorbarticks = isolevels(ctx, funcs)

    if crange[1] == crange[2]
        eps = 1.0e-5
    else
        eps = (crange[2] - crange[1]) * 1.0e-15
    end

    colorlevels = range(crange[1] - eps, crange[2] + eps; length = ctx[:colorlevels])

    planes = makeplanes(xyzmin, xyzmax, ctx[:xplanes], ctx[:yplanes], ctx[:zplanes])

    ccoord0, faces0, values = marching_tetrahedra(
        grids, funcs, planes, levels; gridscale = ctx[:gridscale],
        tol = ctx[:tetxplane_tol]
    )

    faces = reshape(reinterpret(Int32, faces0), (3, length(faces0)))
    ccoord = reshape(reinterpret(Float32, ccoord0), (3, length(ccoord0)))

    nfaces = size(faces, 2)
    if nfaces > 0
        colors = zeros(nfaces)
        for i in 1:nfaces
            colors[i] = (values[faces[1, i]] + values[faces[2, i]] + values[faces[3, i]]) / 3
        end
        # thx, https://stackoverflow.com/a/24229480/8922290
        collec = ctx[:ax].plot_trisurf(
            ccoord[1, :],
            ccoord[2, :],
            transpose(faces .- 1),
            ccoord[3, :];
            cmap = PyPlot.ColorMap(plaincolormap(ctx)),
            vmin = crange[1],
            vmax = crange[2],
        )
        collec.set_array(colors)
        collec.autoscale()
    end

    ax.set_xlim3d(xyzmin[1], xyzmax[1])
    ax.set_ylim3d(xyzmin[2], xyzmax[2])
    ax.set_zlim3d(xyzmin[3], xyzmax[3])
    ax.view_init(ctx[:elev], ctx[:azim])

    if ctx[:colorbar] == :horizontal
        ctx[:cbar] = fig.colorbar(
            collec;
            ax = ax,
            ticks = colorbarticks,
            boundaries = colorlevels,
            orientation = "horizontal",
        )
    end
    if ctx[:colorbar] == :vertical
        ctx[:cbar] = fig.colorbar(
            collec;
            ax = ax,
            ticks = colorbarticks,
            boundaries = colorlevels,
            orientation = "vertical",
        )
    end
    ax.set_title(ctx[:title])

    if ctx[:legend] != :none
        ax.legend(; loc = leglocs[ctx[:legend]])
    end
    return reveal(ctx, TP)
end

### 2D Vector
function vectorplot!(ctx, TP::Type{PyPlotType}, ::Type{Val{2}}, grid, func)
    PyPlot = ctx[:Plotter]
    if !haskey(ctx, :ax)
        ctx[:ax] = ctx[:figure].add_subplot(ctx[:layout]..., ctx[:iplot])
    end
    if ctx[:clear]
        if haskey(ctx, :cbar)
            ctx[:cbar].remove()
        end
        ctx[:ax].remove()
        ctx[:ax] = ctx[:figure].add_subplot(ctx[:layout]..., ctx[:iplot])

        xlimits = ctx[:xlimits]
        ylimits = ctx[:ylimits]
        if xlimits[1] < xlimits[2]
            ctx[:ax].set_xlim(xlimits...)
        end
        if ylimits[1] < ylimits[2]
            ctx[:ax].set_ylim(ylimits...)
        end
    end

    ax = ctx[:ax]
    fig = ctx[:figure]
    ax.set_aspect(ctx[:aspect])
    ax.set_title(ctx[:title])

    rc, rv = vectorsample(
        grid, func; rasterpoints = ctx[:rasterpoints], offset = ctx[:offset], xlimits = ctx[:xlimits],
        ylimits = ctx[:ylimits], gridscale = ctx[:gridscale]
    )
    qc, qv = quiverdata(rc, rv; vscale = ctx[:vscale], vnormalize = ctx[:vnormalize], vconstant = ctx[:vconstant])

    # For the kwargs, see
    # https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.quiver.html
    # Without them, PyPlot itself normalizes
    ax.quiver(qc[1, :], qc[2, :], qv[1, :], qv[2, :]; color = rgbtuple(ctx[:color]), scale = 1, angles = "xy", scale_units = "xy")
    ax.set_xlabel(ctx[:xlabel])
    ax.set_ylabel(ctx[:ylabel])

    return reveal(ctx, TP)
end

### 2D stream
function streamplot!(ctx, TP::Type{PyPlotType}, ::Type{Val{2}}, grid, func)
    PyPlot = ctx[:Plotter]
    if !haskey(ctx, :ax)
        ctx[:ax] = ctx[:figure].add_subplot(ctx[:layout]..., ctx[:iplot])
    end
    if ctx[:clear]
        if haskey(ctx, :cbar)
            ctx[:cbar].remove()
        end
        ctx[:ax].remove()
        ctx[:ax] = ctx[:figure].add_subplot(ctx[:layout]..., ctx[:iplot])
    end

    ax = ctx[:ax]
    fig = ctx[:figure]
    ax.set_aspect(ctx[:aspect])
    ax.set_title(ctx[:title])

    xlimits = ctx[:xlimits]
    ylimits = ctx[:ylimits]
    if xlimits[1] < xlimits[2]
        ax.set_xlim(xlimits...)
    end
    if ylimits[1] < ylimits[2]
        ax.set_ylim(ylimits...)
    end

    # thx https://discourse.julialang.org/t/meshgrid-function-in-julia/48679/4?u=j-fu
    function meshgrid(rc)
        nx = length(rc[1])
        ny = length(rc[2])
        xout = zeros(ny, nx)
        yout = zeros(ny, nx)
        for jx in 1:nx
            for ix in 1:ny
                xout[ix, jx] = rc[1][jx]
                yout[ix, jx] = rc[2][ix]
            end
        end
        return xout, yout
    end

    rc, rv = vectorsample(
        grid, func; rasterpoints = 2 * ctx[:rasterpoints], offset = ctx[:offset], xlimits = ctx[:xlimits],
        ylimits = ctx[:ylimits], gridscale = ctx[:gridscale]
    )

    X, Y = meshgrid(rc)

    ax.streamplot(X, Y, rv[1, :, :, 1]', rv[2, :, :, 1]'; color = rgbtuple(ctx[:color]), density = 1)
    ax.set_xlabel(ctx[:xlabel])
    ax.set_ylabel(ctx[:ylabel])

    return reveal(ctx, TP)
end

function customplot!(ctx, TP::Type{PyPlotType}, func)
    PyPlot = ctx[:Plotter]
    if !haskey(ctx, :ax)
        ctx[:ax] = ctx[:figure].add_subplot(ctx[:layout]..., ctx[:iplot])
    end
    # if ctx[:clear]
    #     if haskey(ctx, :cbar)
    #         ctx[:cbar].remove()
    #     end
    #     ctx[:ax].remove()
    #     ctx[:ax] = ctx[:figure].add_subplot(ctx[:layout]..., ctx[:iplot])
    # end
    func(ctx[:ax])
    return reveal(ctx, TP)
end
