"""
    module GridVisualizeMeshCatExt
    
Extension module for MeshCat.jl. Experimental.
"""
module GridVisualizeMeshCatExt

using Colors
using ColorSchemes
using DocStringExtensions
using GeometryBasics
import GridVisualize: initialize!, save, reveal, gridplot!, scalarplot!, vectorplot!, streamplot!, customplot!
using GridVisualize: MeshCatType, GridVisualizer, SubVisualizer
using GridVisualize: isolevels, cellcolors, num_cellcolors, regionmesh, bfacesegments3
using ExtendableGrids
using GridVisualizeTools

function initialize!(p::GridVisualizer, ::Type{MeshCatType})
    MeshCat = p.context[:Plotter]
    layout = p.context[:layout]
    @assert(layout == (1, 1))
    vis = MeshCat.Visualizer()
    MeshCat.send(vis.core, MeshCat.SetProperty(MeshCat.Path(["Grid"]), "visible", false))
    MeshCat.send(
        vis.core,
        MeshCat.SetProperty(MeshCat.Path(["Background"]), "visible", false)
    )
    p.context[:scene] = vis
    for I in CartesianIndices(layout)
        ctx = p.subplots[I]
        ctx[:figure] = p.context[:scene]
    end
    return nothing
end

function reveal(p::GridVisualizer, ::Type{MeshCatType})
    MeshCat = p.context[:Plotter]
    return MeshCat.IJuliaCell(p.context[:scene])
end

function reveal(ctx::SubVisualizer, TP::Type{MeshCatType})
    if ctx[:show] || ctx[:reveal]
        return reveal(ctx[:GridVisualizer], TP)
    end
    return nothing
end

gridplot!(ctx, TP::Type{MeshCatType}, ::Type{Val{1}}, grid) = nothing
scalarplot!(ctx, TP::Type{MeshCatType}, ::Type{Val{1}}, grids, parentgrid, funcs) = nothing

# 2D grid
function gridplot!(ctx, TP::Type{MeshCatType}, ::Type{Val{2}}, grid)
    MeshCat = ctx[:Plotter]
    vis = ctx[:figure]

    nregions = num_cellregions(grid)
    nbregions = num_bfaceregions(grid)

    cmap = region_cmap(nregions)
    bcmap = bregion_cmap(nbregions)
    for i in 1:nregions
        mesh = regionmesh(grid, 1.0, i; cellcoloring = ctx[:cellcoloring])
        MeshCat.setobject!(
            vis["interior"]["r$(i)"],
            mesh,
            MeshCat.MeshLambertMaterial(; color = RGBA{Float32}(cmap[i], 1.0))
        )
        MeshCat.setobject!(
            vis["interior"]["r$(i)_edges"],
            mesh,
            MeshCat.MeshPhongMaterial(;
                color = RGBA{Float32}(0.0, 0.0, 0.0, 1.0),
                wireframe = true,
            )
        )
    end

    for i in 1:nbregions
        points = bfacesegments3(grid, 1.0, i)
        mat = MeshCat.MeshLambertMaterial(; color = RGBA{Float32}(bcmap[i], 1.0))
        ls = MeshCat.LineSegments(points, mat)
        MeshCat.setobject!(vis["boundary"]["b$(i)"], ls)
    end
    MeshCat.send(vis.core, MeshCat.SetProperty(MeshCat.Path(["Axes"]), "visible", false))

    return reveal(ctx, TP)
end

function gridplot!(ctx, TP::Type{MeshCatType}, ::Type{Val{3}}, grid)
    MeshCat = ctx[:Plotter]
    vis = ctx[:figure]

    nregions = num_cellregions(grid)
    nbregions = num_bfaceregions(grid)
    cmap = region_cmap(nregions)
    bcmap = bregion_cmap(nbregions)

    xyzmin = zeros(3)
    xyzmax = ones(3)
    coord = grid[Coordinates]
    @views for idim in 1:3
        xyzmin[idim] = minimum(coord[idim, :])
        xyzmax[idim] = maximum(coord[idim, :])
    end

    ctx[:xplane] = max(xyzmin[1], min(xyzmax[1], ctx[:xplanes][1]))
    ctx[:yplane] = max(xyzmin[2], min(xyzmax[2], ctx[:yplanes][1]))
    ctx[:zplane] = max(xyzmin[3], min(xyzmax[3], ctx[:zplanes][1]))

    xyzcut = [ctx[:xplane], ctx[:yplane], ctx[:zplane]]

    if ctx[:interior]
        pts, fcs = extract_visible_cells3D(
            grid,
            xyzcut;
            cellcoloring = ctx[:cellcoloring],
            primepoints = hcat(xyzmin, xyzmax),
            Tp = Point3f,
            Tf = GLTriangleFace,
        )

        for i in 1:nregions
            mesh = Mesh(pts[i], fcs[i])
            MeshCat.setobject!(
                vis["r$(i)"],
                mesh,
                MeshCat.MeshLambertMaterial(; color = RGBA{Float32}(cmap[i], 1.0))
            )
            MeshCat.setobject!(
                vis["r$(i)_edges"],
                mesh,
                MeshCat.MeshPhongMaterial(;
                    color = RGBA{Float32}(0.0, 0.0, 0.0, 1.0),
                    wireframe = true,
                )
            )
        end
    end

    pts, fcs = extract_visible_bfaces3D(
        grid,
        xyzcut;
        primepoints = hcat(xyzmin, xyzmax),
        Tp = Point3f,
        Tf = GLTriangleFace,
    )

    for i in 1:nbregions
        mesh = Mesh(pts[i], fcs[i])
        MeshCat.setobject!(
            vis["b$(i)"],
            mesh,
            MeshCat.MeshLambertMaterial(; color = RGBA{Float32}(bcmap[i], 1.0))
        )
        MeshCat.setobject!(
            vis["b$(i)_edges"],
            mesh,
            MeshCat.MeshPhongMaterial(;
                color = RGBA{Float32}(0.0, 0.0, 0.0, 1.0),
                wireframe = true,
            )
        )
    end

    return reveal(ctx, TP)
end

function scalarplot!(ctx, TP::Type{MeshCatType}, ::Type{Val{3}}, grids, parentgrid, funcs)
    grid = parentgrid
    func = funcs[1]
    MeshCat = ctx[:Plotter]
    vis = ctx[:figure]

    nregions = num_cellregions(grid)
    nbregions = num_bfaceregions(grid)
    bcmap = bregion_cmap(nbregions)
    xyzmin = zeros(3)
    xyzmax = ones(3)
    coord = grid[Coordinates]
    @views for idim in 1:3
        xyzmin[idim] = minimum(coord[idim, :])
        xyzmax[idim] = maximum(coord[idim, :])
    end
    xyzcut = [ctx[:xplanes][1], ctx[:yplanes][1], ctx[:zplanes][1]]
    fminmax = extrema(func)

    ctx[:xplane] = max(xyzmin[1], min(xyzmax[1], ctx[:xplanes][1]))
    ctx[:yplane] = max(xyzmin[2], min(xyzmax[2], ctx[:yplanes][1]))
    ctx[:zplane] = max(xyzmin[3], min(xyzmax[3], ctx[:zplanes][1]))
    ctx[:flevel] = max(fminmax[1], min(fminmax[2], ctx[:levels][1]))

    makeplanes(x, y, z) = [[1, 0, 0, -x], [0, 1, 0, -y], [0, 0, 1, -z]]

    ccoord, faces, values = marching_tetrahedra(
        grid,
        func,
        makeplanes(ctx[:xplane], ctx[:yplane], ctx[:zplane]),
        [ctx[:flevel]];
        primepoints = hcat(xyzmin, xyzmax),
        primevalues = fminmax,
        Tp = Point3f,
        Tf = GLTriangleFace,
        Tv = Float32,
    )
    mesh = Mesh(ccoord, faces)

    to01(v) = (v - fminmax[1]) / (fminmax[2] - fminmax[1])
    rgb(v) = RGB(v, 0.0, 1.0 - v)

    vcmap = colorschemes[ctx[:colormap]]
    mesh = GeometryBasics.Mesh(ccoord, faces; color = [get(vcmap, values[i]) for i in 1:length(values)], normal = normals(ccoord, faces))
    material = MeshCat.MeshLambertMaterial(; vertexColors = true)

    MeshCat.setobject!(vis[:marching_tets], mesh, material)

    if ctx[:outlinealpha] > 0
        pts, fcs = extract_visible_bfaces3D(
            grid,
            xyzmax;
            primepoints = hcat(xyzmin, xyzmax),
            Tp = Point3f,
            Tf = GLTriangleFace,
        )
        for i in 1:nbregions
            mesh = Mesh(pts[i], fcs[i])
            MeshCat.setobject!(
                vis["b$(i)"],
                mesh,
                MeshCat.MeshLambertMaterial(;
                    color = color = RGBA{Float32}(bcmap[i], 0.35),
                )
            )
        end
    end
    return reveal(ctx, TP)
end
end
