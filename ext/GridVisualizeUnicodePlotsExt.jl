"""
    module GridVisualizeUnicodePlotsExt

Extension module for UnicodePlots.jl
"""
module GridVisualizeUnicodePlotsExt

import GridVisualize: initialize!, gridplot!, scalarplot!, vectorplot!, bregion_cmap, region_cmap, reveal, streamplot!, compute_3d_z_rotation_matrix
using GridVisualize: UnicodePlotsType, GridVisualizer, SubVisualizer, vectorsample, quiverdata, isolevels
using GridVisualizeTools: marching_triangles
using UnicodePlots: UnicodePlots
using ExtendableGrids: Coordinates, simplexgrid, ON_CELLS, ON_FACES, ON_EDGES, CellNodes, FaceNodes, BFaceNodes, CellGeometries, CellRegions, BFaceRegions, num_cells, num_nodes, local_celledgenodes, num_bfaceregions, num_cellregions, num_targets, interpolate!
using Colors: Colors, RGB, RGBA, red, green, blue
using ColorSchemes: colorschemes, color

initialize!(p, ::Type{UnicodePlotsType}) = nothing


function reveal(p::GridVisualizer, TP::Type{UnicodePlotsType})
    layout = p.context[:layout]
    subplots = @views permutedims(p.subplots)[:]

    if layout == (1, 1)
        display(subplots[1][:figure])
    else
        if :Term ∉ Symbol.(Base.loaded_modules_array())
            @warn "A GridVisualizer with multiple UnicodePlots can only respect the desired layout if 'Term.jl' is loaded: add Term.jl to your environment and use it."

            # plot figures one by one
            for subplot in subplots
                haskey(subplot, :figure) && display(subplot[:figure])
            end
        else
            # set missing figures
            for subplot in subplots
                if !haskey(subplot, :figure)
                    missing_figure!(subplot, TP)
                end
            end

            grid_plot = UnicodePlots.gridplot([subplot[:figure] for subplot in subplots], layout = p.context[:layout], show_placeholder = true)
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

# unit vector from origin towards the camera for given azimuth/elevation angles (degrees)
function camera_direction(azim::Real, elev::Real)
    az = Float64(azim) * π / 180
    el = Float64(elev) * π / 180
    return [cos(el) * cos(az), cos(el) * sin(az), sin(el)]
end

# cyclic polygon edges of face j (as node index pairs)
function face_edges(facenodes, j)
    n = num_targets(facenodes, j)
    return [(facenodes[k, j], facenodes[k % n + 1, j]) for k in 1:n]
end

function gridplot!(ctx, TP::Type{UnicodePlotsType}, ::Type{Val{3}}, grid)
    UnicodePlots = ctx[:Plotter]

    # set up camera from azimuth / elevation:
    # R rotates the z axis into the view direction (pointing towards the camera),
    # so in rotated coordinates q = R'p the first two components are screen
    # coordinates and the third one is the depth (larger = closer to camera).
    # Coordinates are centered w.r.t. the bounding box center.
    R = compute_3d_z_rotation_matrix(camera_direction(ctx[:azim], ctx[:elev]))
    coords = grid[Coordinates]
    c3 = vec((minimum(coords, dims = 2) .+ maximum(coords, dims = 2)) / 2)
    kappa = Float64(ctx[:perspectiveness])

    q = R' * (coords .- c3)
    d0, d1 = extrema(q[3, :])
    ddp = max(d1 - d0, 1.0e-14)

    # perspective: distant points shrink towards the view center
    if kappa > 0
        sc = reshape(1 ./ (1 .+ kappa .* clamp.((d1 .- q[3:3, :]) ./ ddp, 0.0, 1.0)), 1, :)
        pc = [q[1:1, :] .* sc; q[2:2, :] .* sc; q[3:3, :]]
    else
        pc = q
    end
    eu = extrema(view(pc, 1, :))
    ev = extrema(view(pc, 2, :))

    # project an arbitrary 3D point (same transform, including perspective)
    function project3(p)
        r = R' * (p .- c3)
        if kappa > 0
            s = 1 / (1 + kappa * clamp((d1 - r[3]) / ddp, 0.0, 1.0))
            return [r[1] * s, r[2] * s, r[3]]
        end
        return r
    end

    # determine resolution (divided by 10, to reduce pixel count in the terminal)
    layout = ctx[:layout]
    resolution = ctx[:size] ./ (48, 24) ./ (layout[2], layout[1])
    aspect = ctx[:aspect]

    # ensure that legend fits
    ncellregions = num_cellregions(grid)
    nbregions = num_bfaceregions(grid)
    resolution = (resolution[1], max(resolution[2], 5 + ncellregions + nbregions))

    # rescale resolution
    wu = eu[2] - eu[1]
    wv = ev[2] - ev[1]
    rescale = 2 * wu / wv * (resolution[2] / (resolution[1])) / aspect
    resolution = (resolution[1] * rescale, resolution[2])

    # we need an integer resolution
    resolution = @. Int(round(resolution))

    # reserve an empty band below the mesh s.t. the axis triad sits outside the geometry
    # (UnicodePlots cannot paint outside of the plot frame, so an empty band is used)
    nb = max(3, Int(round(0.2 * resolution[2]))) # triad band rows
    vunit = wv / resolution[2]                    # virtual units per row (mesh scale)
    band_h = nb * vunit                           # band height in virtual units
    ev_canvas = (ev[1] - band_h, ev[2])           # extended vertical range
    canvas_rows = resolution[2] + nb

    # create UnicodePlots.Canvas
    CanvasType = UnicodePlots.BrailleCanvas # should this be a changeable parameter ?
    canvas = CanvasType(
        canvas_rows, resolution[1],                   # number of rows and columns (characters)
        origin_y = ev_canvas[1], origin_x = eu[1],    # position in virtual space (incl. band)
        height = (ev_canvas[2] - ev_canvas[1]), width = eu[2] - eu[1]; blend = false
    )

    # region colors
    cmap = region_cmap(max(2, ncellregions))
    ctx[:cmap] = cmap
    cell_colors = [(Int(round(cmap[i].r * 255)), Int(round(cmap[i].g * 255)), Int(round(cmap[i].b * 255))) for i in 1:ncellregions]
    bcmap = bregion_cmap(nbregions)
    ctx[:bcmap] = bcmap
    bcolors = [(Int(round(bcmap[i].r * 255)), Int(round(bcmap[i].g * 255)), Int(round(bcmap[i].b * 255))) for i in 1:nbregions]

    linewidth = ctx[:linewidth]
    if linewidth > 0
        # interior edges (from cells, colored by cell region): painter's algorithm
        # (far edges first), brightness fades with depth
        if ctx[:interior] && ctx[:wireframe]
            cellnodes = grid[CellNodes]
            cellregions = grid[CellRegions]
            cellgeoms = grid[CellGeometries]
            ncells = num_cells(grid)
            ed1 = Int[]
            ed2 = Int[]
            edd = Float64[]
            edr = Int[]
            for j in 1:ncells
                cen = local_celledgenodes(cellgeoms[j])
                for k in 1:size(cen, 2)
                    i1 = cellnodes[cen[1, k], j]
                    i2 = cellnodes[cen[2, k], j]
                    push!(ed1, i1)
                    push!(ed2, i2)
                    push!(edd, 0.5 * (pc[3, i1] + pc[3, i2]))
                    push!(edr, cellregions[j])
                end
            end
            if !isempty(edd)
                dmin, dmax = extrema(edd)
                dd = dmax - dmin
                dd == 0 && (dd = 1.0)
                for k in sortperm(edd)
                    # brightness encodes depth: closer edges are brighter
                    s = 0.35 + 0.65 * (edd[k] - dmin) / dd
                    c0 = cell_colors[edr[k]]
                    edge_color = (Int(round(c0[1] * s)), Int(round(c0[2] * s)), Int(round(c0[3] * s)))
                    UnicodePlots.lines!(
                        canvas,
                        pc[1, ed1[k]], pc[2, ed1[k]], # from
                        pc[1, ed2[k]], pc[2, ed2[k]]; # to
                        color = edge_color
                    )
                end
            end
        end

        # boundary faces on top, colored by bregion, with depth fog
        bfacenodes = grid[BFaceNodes]
        bfaceregions = grid[BFaceRegions]
        nbfaces = size(bfacenodes, 2)
        bdepth = [sum(pc[3, bfacenodes[1:num_targets(bfacenodes, j), j]]) / num_targets(bfacenodes, j) for j in 1:nbfaces]
        bd0, bd1 = extrema(bdepth)
        bdd = max(bd1 - bd0, 1.0e-14)
        for j in sortperm(bdepth)
            s = 0.45 + 0.55 * (bdepth[j] - bd0) / bdd
            c0 = bcolors[bfaceregions[j]]
            edge_color = (Int(round(c0[1] * s)), Int(round(c0[2] * s)), Int(round(c0[3] * s)))
            for (i1, i2) in face_edges(bfacenodes, j)
                UnicodePlots.lines!(
                    canvas,
                    pc[1, i1], pc[2, i1], # from
                    pc[1, i2], pc[2, i2]; # to
                    color = edge_color
                )
            end
        end
    end

    # mark cell midpoints with their region number in region color (opt-in via cellmarkers)
    markersize = ctx[:markersize]
    if ctx[:interior] && ctx[:wireframe] && ctx[:cellmarkers]
        cellregions = grid[CellRegions]
        cellnodes = grid[CellNodes]
        ncells = num_cells(grid)
        cells = Tuple{Float64, Float64, Float64, Int}[] # (depth, u, v, region)
        midpoint = zeros(3)
        p2 = zeros(3)
        for j in 1:ncells
            fill!(midpoint, 0.0)
            nvertices = num_targets(cellnodes, j)
            for k in 1:nvertices
                midpoint .+= coords[:, cellnodes[k, j]]
            end
            midpoint ./= nvertices
            p2 .= project3(midpoint)
            push!(cells, (p2[3], p2[1], p2[2], cellregions[j]))
        end
        # painter's algorithm: far cells first, closer numbers overwrite them
        sort!(cells; by = c -> c[1])
        for c in cells
            r = c[4]
            if markersize < 4
                UnicodePlots.points!(
                    canvas,
                    c[2], c[3];
                    color = cell_colors[r]
                )
            else
                UnicodePlots.annotate!(
                    canvas,
                    c[2], c[3],
                    string(r),
                    UnicodePlots.ansi_color(cell_colors[r]),
                    false
                )
            end
        end
    end

    plt = UnicodePlots.Plot(canvas; title = ctx[:title], border = ctx[:border])

    # coordinate hints for the projected view:
    # (a) world coordinates at the bounding box corner regions, snapped to the closest
    # mesh node s.t. labels and markers sit on actual geometry (near labels win)
    lo = vec(minimum(coords, dims = 2))
    hi = vec(maximum(coords, dims = 2))
    corners = [[x, y, z] for x in (lo[1], hi[1]) for y in (lo[2], hi[2]) for z in (lo[3], hi[3])]
    fmt(x) = isinteger(x) ? string(Int(x)) : string(round(x, sigdigits = 3))
    dcam = R[:, 3] # view direction (pointing towards the camera)
    nnodes = num_nodes(grid)
    hints = Tuple{Float64, Float64, Float64, Int}[] # (depth, u, v, node)
    seen_nodes = Int[]
    for c in corners
        # consider only corner regions which are really visible, i.e. at least one of
        # the three adjacent bounding box faces points towards the camera
        any((c[i] == hi[i] ? 1 : -1) * dcam[i] > 0 for i in 1:3) || continue
        # snap to the closest mesh node, skip corners sharing the same node
        jnode = 1
        dmin = Inf
        for i in 1:nnodes
            d = (coords[1, i] - c[1])^2 + (coords[2, i] - c[2])^2 + (coords[3, i] - c[3])^2
            d < dmin && (dmin = d; jnode = i)
        end
        jnode in seen_nodes && continue
        push!(seen_nodes, jnode)
        pp = project3(view(coords, :, jnode))
        push!(hints, (pp[3], pp[1], pp[2], jnode))
    end
    placed = Tuple{Float64, Float64}[]
    for h in sort(hints; by = h -> h[1]) # far nodes first, near labels overwrite them
        (u, v, jnode) = (h[2], h[3], h[4])
        pn = view(coords, :, jnode)
        # skip labels colliding with an already placed one (approximately 9 chars wide)
        if !any(abs(u - p0[1]) < 9 * wu / resolution[1] && abs(v - p0[2]) < 2 * wv / resolution[2] for p0 in placed)
            push!(placed, (u, v))
            text = "($(fmt(pn[1])), $(fmt(pn[2])), $(fmt(pn[3])))"
            # shift towards the canvas center s.t. labels do not get clipped at the border
            half = (length(text) / 2 + 1) * wu / resolution[1]
            du = u > (eu[1] + eu[2]) / 2 ? -half : half
            dv = v > (ev[1] + ev[2]) / 2 ? -1.5 * wv / resolution[2] : 0.5 * wv / resolution[2]
            # bullet marker at the annotated node itself, slightly shifted inwards
            UnicodePlots.annotate!(
                canvas,
                u + (du > 0 ? 0.7 : -0.7) * wu / resolution[1],
                v + (dv > 0 ? 0.7 : -0.7) * wv / resolution[2],
                "●",
                UnicodePlots.ansi_color((230, 230, 230)),
                false
            )
            UnicodePlots.annotate!(
                canvas, u + du, v + dv,
                text,
                UnicodePlots.ansi_color((200, 200, 200)),
                false
            )
        end
    end

    # (b) projected axis triad, neutral gray (to not be confused with region colors),
    # placed in the empty band reserved below the mesh, drawn last
    L = min(0.42 * band_h, 0.4 * wu)
    tips = [(R' * (L * [Float64(k) for k in (i, j, l)]))[1:2] for (i, j, l) in ((1, 0, 0), (0, 1, 0), (0, 0, 1))]
    # anchor the triad in the bottom-left of the reserved band, all arrows inside
    lo_v = min(0.0, minimum(t[2] for t in tips))
    hi_v = max(0.0, maximum(t[2] for t in tips))
    tx = eu[1] + 0.06 * wu - min(0.0, minimum(t[1] for t in tips))
    ty = ev_canvas[1] + (band_h - (hi_v - lo_v)) / 2 - lo_v
    triad_color = (180, 180, 180)
    for (t2, axis_label) in zip(tips, "xyz")
        UnicodePlots.lines!(canvas, tx, ty, tx + t2[1], ty + t2[2]; color = triad_color)
        # place letter at the arrow tip, slightly offset to the right s.t. it is not clipped
        UnicodePlots.annotate!(
            canvas,
            tx + t2[1] + 0.8 * wu / resolution[1], ty + t2[2] + 0.5 * wv / resolution[2],
            string(axis_label),
            UnicodePlots.ansi_color(triad_color),
            false
        )
    end

    y0 = 0
    if ctx[:interior] && ctx[:wireframe] && ctx[:cellmarkers]
        y0 = region_legend!(plt, " cell", 1, [])
        y0 = region_legend!(plt, "regions", 2, cell_colors)
    end
    region_legend!(plt, " bface", y0 + 2, [])
    region_legend!(plt, "regions", y0 + 3, bcolors)

    # camera info (view space u/v extents would be meaningless to the user,
    # world coordinates are annotated at the corner nodes instead)
    UnicodePlots.label!(plt, :b, "azim = $(ctx[:azim])°, elev = $(ctx[:elev])°")

    # plot
    ctx[:figure] = plt
    return reveal(ctx, TP)
end

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
    resolution = ctx[:size] ./ (48, 24) ./ (layout[2], layout[1])
    aspect = ctx[:aspect]

    # ensure that legend fits
    ncellregions = num_cellregions(grid)
    nbregions = num_bfaceregions(grid)
    resolution = (resolution[1], max(resolution[2], 5 + ncellregions + nbregions))

    # rescale resolution
    wx = ex[2] - ex[1]
    wy = ey[2] - ey[1]
    rescale = 2 * wx / wy * (resolution[2] / (resolution[1])) / aspect
    resolution = (resolution[1] * rescale, resolution[2])

    # we need an integer resolution
    resolution = @. Int(round(resolution))

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
    cell_colors = [(Int(round(cmap[i].r * 255)), Int(round(cmap[i].g * 255)), Int(round(cmap[i].b * 255))) for i in 1:ncellregions]
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
    bcolors = [(Int(round(bcmap[i].r * 255)), Int(round(bcmap[i].g * 255)), Int(round(bcmap[i].b * 255))) for i in 1:nbregions]
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
    resolution = (Int(round(ctx[:size][1] / 12 / layout[2])), max(7, 5 + ncellregions + nbregions))

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
    colors = [(Int(round(cmap[i].r * 255)), Int(round(cmap[i].g * 255)), Int(round(cmap[i].b * 255))) for i in 1:ncellregions]
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
    bcolors = [(Int(round(bcmap[i].r * 255)), Int(round(bcmap[i].g * 255)), Int(round(bcmap[i].b * 255))) for i in 1:nbregions]
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
    resolution = @. Int(round(ctx[:size] ./ 12 ./ (layout[2], 2 * layout[1]))) # reduce pixel count in the terminal (size is then compatible to other plots)

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
    resolution = ctx[:size] ./ 12 ./ (layout[2], layout[1]) # reduce pixel count in the terminal
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

    # rescale resolution
    wx = ex[2] - ex[1]
    wy = ey[2] - ey[1]
    rescale = wx / wy * (resolution[2] / (resolution[1])) / aspect
    resolution = (resolution[1] * rescale, resolution[2])

    # we need an integer resolution
    resolution = @. Int(round(resolution))

    X = LinRange(ex[1], ex[2], resolution[1])
    Y = LinRange(ey[1], ey[2], resolution[2])
    xgrid_plot = simplexgrid(X, Y)

    # interpolate data onto plot_grid
    I = zeros(Float64, num_nodes(xgrid_plot))
    interpolate!(I, xgrid_plot, func, grids[1]; eps = 1.0e-14, not_in_domain_value = NaN, trybrute = true)

    @info extrema(func) size (func)

    # adjust colormap to desired number of colorlevels
    colormap = get(colorschemes[ctx[:colormap]], range(0.0, 1.0, length = max(2, 1 + ctx[:colorlevels])))

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

    # isolines (only when given as an array)
    if typeof(ctx[:levels]) <: AbstractArray
        levels, ~, ~ = isolevels(ctx, funcs)
        points = marching_triangles(grids[1], func, levels; gridscale = 1.0)
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
        for j in 1:2:(length(points) - 1)
            UnicodePlots.lineplot!(plt, [points[j][1], points[j + 1][1]], [points[j][2], points[j + 1][2]], color = color)
        end
    end

    ctx[:figure] = plt

    return reveal(ctx, TP)
end

scalarplot!(ctx, TP::Type{UnicodePlotsType}, ::Type{Val{3}}, grids, parentgrid, funcs) = missing_figure!(ctx, TP, "scalarplot! (3D)", not_implemented = true)


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
    resolution = ctx[:size] ./ (48, 24) ./ (layout[2], layout[1])

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

    # rescale resolution
    wx = ex[2] - ex[1]
    wy = ey[2] - ey[1]
    rescale = 2 * wx / wy * (resolution[2] / (resolution[1])) / aspect
    resolution = (resolution[1] * rescale, resolution[2])

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

streamplot!(ctx, TP::Type{UnicodePlotsType}, ::Type{Val{2}}, grid, func) = missing_figure!(ctx, TP, "streamplot! (2D)", not_implemented = true)

# placeholder for missing/not implemented
function missing_figure!(ctx, TP::Type{UnicodePlotsType}, caller_text = ""; not_implemented = false)

    UnicodePlots = ctx[:Plotter]

    layout = ctx[:layout]
    resolution = ctx[:size] ./ (12, 24) ./ (layout[2], layout[1])

    # we need an integer resolution
    resolution = @. Int(round(resolution))

    # create UnicodePlots.Canvas
    CanvasType = UnicodePlots.BrailleCanvas # should this be a changeable parameter ?
    canvas = CanvasType(
        resolution[2], resolution[1]            # number of rows and columns (characters)
    )

    if not_implemented
        @warn "$caller_text is not implemented for the UnicodePlots backend"
        UnicodePlots.annotate!(
            canvas,
            0.5, 0.4,
            "[not implemented]",
            UnicodePlots.ansi_color(:normal),
            false
        )
    end

    UnicodePlots.annotate!(
        canvas,
        0.5,
        not_implemented ? 0.6 : 0.5,
        caller_text,
        UnicodePlots.ansi_color(:normal),
        false
    )

    plt = UnicodePlots.Plot(canvas; title = ctx[:title], border = :dotted)

    # plot
    ctx[:figure] = plt
    return reveal(ctx, TP)
end

end # module
