"""
$(SIGNATURES)

Extract visible tetrahedra - those intersecting with the planes
`x=xyzcut[1]` or `y=xyzcut[2]`  or `z=xyzcut[3]`. 

Return corresponding points and facets for each region for drawing as mesh (Makie,MeshCat)
or trisurf (pyplot)
"""
function GridVisualizeTools.extract_visible_cells3D(
        grid::ExtendableGrid, xyzcut;
        cellcoloring = :cellregions,
        gridscale = 1.0,
        primepoints = zeros(0, 0),
        Tp = SVector{3, Float32},
        Tf = SVector{3, Int32}
    )
    coord = grid[Coordinates] * gridscale
    cellnodes = grid[CellNodes]
    cellregions = cellcolors(grid, cellcoloring)
    nregions = num_cellcolors(grid, cellcoloring)
    return extract_visible_cells3D(
        coord, cellnodes, cellregions, nregions, [xyzcut...] * gridscale;
        primepoints = primepoints,
        Tp = Tp, Tf = Tf
    )
end

"""
$(SIGNATURES)

Extract visible boundary faces - those not cut off by the planes
`x=xyzcut[1]` or `y=xyzcut[2]`  or `z=xyzcut[3]`. 

Return corresponding points and facets for each region for drawing as mesh (Makie,MeshCat)
or trisurf (pyplot)
"""
function GridVisualizeTools.extract_visible_bfaces3D(
        grid::ExtendableGrid, xyzcut;
        gridscale = 1.0,
        primepoints = zeros(0, 0),
        Tp = SVector{3, Float32},
        Tf = SVector{3, Int32}
    )
    coord = grid[Coordinates] * gridscale
    bfacenodes = grid[BFaceNodes]
    bfaceregions = grid[BFaceRegions]
    nbregions = grid[NumBFaceRegions]

    return extract_visible_bfaces3D(
        coord, bfacenodes, bfaceregions, nbregions, [xyzcut...] * gridscale;
        primepoints = primepoints, Tp = Tp, Tf = Tf
    )
end

"""
   $(SIGNATURES)

Extract isosurfaces and plane interpolation for function on 3D tetrahedral mesh.
See [`marching_tetrahedra(coord,cellnodes,func,planes,flevels;tol, primepoints, primevalues, Tv, Tp, Tf)`](@ref)
"""
function GridVisualizeTools.marching_tetrahedra(
        grid::ExtendableGrid, func, planes, flevels; gridscale = 1.0,
        kwargs...
    )
    coord = grid[Coordinates] * gridscale
    cellnodes = grid[CellNodes]
    return marching_tetrahedra(coord, cellnodes, func, planes, flevels; kwargs...)
end

function GridVisualizeTools.marching_tetrahedra(
        grids::Vector{ExtendableGrid{Tv, Ti}}, funcs, planes, flevels;
        gridscale = 1.0,
        kwargs...
    ) where {Tv, Ti}
    coord = [grid[Coordinates] * gridscale for grid in grids]
    cellnodes = [grid[CellNodes] for grid in grids]
    return marching_tetrahedra(coord, cellnodes, funcs, planes, flevels; kwargs...)
end

########################################################################################
"""
    $(SIGNATURES)

Deprecated
"""
function GridVisualizeTools.marching_triangles(grid::ExtendableGrid, func, levels; gridscale = 1.0)
    coord::Matrix{Float64} = grid[Coordinates] * gridscale
    cellnodes::Matrix{Int32} = grid[CellNodes]
    points, _, _ = marching_triangles(coord, cellnodes, func, [], levels)
    return points
end

"""
    $(SIGNATURES)

Collect isoline snippets and/or intersection points with lines and values ready for linesegments!
"""
function GridVisualizeTools.marching_triangles(grid::ExtendableGrid, func, lines, levels; gridscale = 1.0)
    coord::Matrix{Float64} = grid[Coordinates] * gridscale
    cellnodes::Matrix{Int32} = grid[CellNodes]
    return marching_triangles(coord, cellnodes, func, lines, levels)
end

function GridVisualizeTools.marching_triangles(grids::Vector{ExtendableGrid{Tv, Ti}}, funcs, lines, levels; gridscale = 1.0) where {Tv, Ti}
    coords = [grid[Coordinates] * gridscale for grid in grids]
    cellnodes = [grid[CellNodes] for grid in grids]
    return marching_triangles(coords, cellnodes, funcs, lines, levels)
end

##############################################
# Create meshes from grid data
function regionmesh(grid, gridscale, iregion; cellcoloring = :cellregions)
    coord = grid[Coordinates]
    cn = grid[CellNodes]
    cr = cellcolors(grid, cellcoloring)
    @views points = [Point2f(coord[:, i] * gridscale) for i in 1:size(coord, 2)]
    faces = Vector{GLTriangleFace}(undef, 0)
    for i in 1:length(cr)
        if cr[i] == iregion
            @views push!(faces, cn[:, i])
        end
    end
    return Mesh(points, faces)
end

function bfacesegments(grid, gridscale, ibreg)
    coord = grid[Coordinates]
    nbfaces = num_bfaces(grid)
    bfacenodes = grid[BFaceNodes]
    bfaceregions = grid[BFaceRegions]
    points = Vector{Point2f}(undef, 0)
    for ibface in 1:nbfaces
        if bfaceregions[ibface] == ibreg
            push!(
                points,
                Point2f(coord[1, bfacenodes[1, ibface]] * gridscale, coord[2, bfacenodes[1, ibface]] * gridscale)
            )
            push!(
                points,
                Point2f(coord[1, bfacenodes[2, ibface]] * gridscale, coord[2, bfacenodes[2, ibface]] * gridscale)
            )
        end
    end
    return points
end

function bfacesegments3(grid, gridscale, ibreg)
    coord = grid[Coordinates]
    nbfaces = num_bfaces(grid)
    bfacenodes = grid[BFaceNodes]
    bfaceregions = grid[BFaceRegions]
    points = Vector{Point3f}(undef, 0)
    for ibface in 1:nbfaces
        if bfaceregions[ibface] == ibreg
            push!(
                points,
                Point3f(
                    coord[1, bfacenodes[1, ibface]] * gridscale, coord[2, bfacenodes[1, ibface]] * gridscale,
                    0.0
                )
            )
            push!(
                points,
                Point3f(
                    coord[1, bfacenodes[2, ibface]] * gridscale, coord[2, bfacenodes[2, ibface]] * gridscale,
                    0.0
                )
            )
        end
    end
    return points
end

"""
          vectorsample(grid::ExtendableGrid{Tv, Ti}, v;
                      offset = :default,
                      rasterpoints = 15,
                      reltol = 1.0e-10,
                      gridscale = 1.0,
                      xlimits = (1, -1),
                      ylimits = (1, -1),
                      zlimits = (1, -1)) where {Tv, Ti}

Extract values of given vector field (either nodal values of a piecewise linear vector field or a callback function
providing evaluation of the vector field for given generalized barycentric coordinates).
at all sampling points on  `offset+ i*spacing` for i in Z^d  defined by the tuples offset and spacing.
Returned values can be fed into [`quiverdata`](@ref)

By default, offset is at the minimum of grid coordinates, and spacing is defined
the largest grid extend divided by 10.

The intermediate `rasterflux` in future versions can be used to calculate
streamlines.
    
The code is 3D ready.
"""
function vectorsample(
        grid::ExtendableGrid{Tv, Ti}, v;
        offset = :default,
        rasterpoints = 32,
        reltol = 1.0e-10,
        gridscale = 1.0,
        xlimits = (1, -1),
        ylimits = (1, -1),
        zlimits = (1, -1)
    ) where {Tv, Ti}
    coord = grid[Coordinates] * gridscale
    cn = grid[CellNodes]
    ncells::Int = num_cells(grid)

    dim::Int = dim_space(grid)

    eltype = dim == 2 ? Triangle2D : Tetrahedron3D

    scaledgrid = grid
    if gridscale != 1.0
        scaledgrid = ExtendableGrid{coord_type(grid), index_type(grid)}()
        for (k, v) in grid.components
            scaledgrid.components[k] = v
        end
        scaledgrid[Coordinates] = coord
    end
    L2G = ExtendableGrids.L2GTransformer(eltype, scaledgrid, ON_CELLS)

    # memory for  inverse of local transformation matrix
    invA = zeros(dim + 1, dim + 1)

    # barycentric coordinates
    λ = zeros(dim + 1)

    # coordinate window
    cminmax = extrema(coord; dims = (2,))
    if xlimits[1] < xlimits[2]
        cminmax[1] = xlimits[:]
    end
    if ylimits[1] < ylimits[2]
        cminmax[2] = ylimits[:]
    end
    if zlimits[1] < zlimits[2]
        cminmax[3] = zlimits[:]
    end
    if offset == :default
        offset = [cminmax[i][1] for i in 1:dim]
    else
        offset = offset * gridscale
    end

    # extent  of domain
    extent = maximum([cminmax[i][2] - cminmax[i][1] for i in 1:dim])

    # scale tolerance
    tol = reltol * extent

    # point spacing
    spacing = [extent / rasterpoints for i in 1:dim]

    # index range
    ijkmax = ones(Int, 3)
    for idim in 1:dim
        ijkmax[idim] = ceil(Int64, (cminmax[idim][2] - offset[idim]) / spacing[idim]) + 1
    end

    # The ijk raster corresponds to a  tensorproduct grid
    # spanned by x,y and z coordinate vectors. Here, we build them
    # in order to avoid to calculate them from the raster indices
    rastercoord = [zeros(Float32, ijkmax[idim]) for idim in 1:dim]

    for idim in 1:dim
        rastercoord[idim] = collect(
            range(
                offset[idim]; step = spacing[idim],
                length = ijkmax[idim]
            )
        )
    end

    # Memory for flux vectors on ijk grid
    rasterflux = zeros(Float32, dim, ijkmax[1], ijkmax[2], ijkmax[3])

    # type stable versions of offset and spacing
    O = zeros(dim)
    O .= offset
    S = zeros(dim)
    S .= spacing

    # memory for point, vector to be investigated
    X = zeros(Float32, dim)
    V = zeros(Float32, dim)

    tmin = zeros(Float64, dim)
    tmax = zeros(Float64, dim)

    # index vector
    I = ones(Int, 3)

    # cell extent
    Imin = ones(Int, 3)
    Imax = ones(Int, 3)

    for icell::Int in 1:ncells
        update_trafo!(L2G, icell) # 1 alloc: the only one left in this cell loop

        # Cell coordinate window
        inode = cn[1, icell]
        @views tmin .= coord[:, inode]
        @views tmax .= coord[:, inode]
        for i in 2:(dim + 1)
            inode = cn[i, icell]
            for idim in 1:dim
                tmin[idim] = min(tmin[idim], coord[idim, inode])
                tmax[idim] = max(tmax[idim], coord[idim, inode])
            end
        end

        # min and max of raster indices falling into cell coordinate window
        for idim in 1:dim
            Imin[idim] = floor(Int64, (tmin[idim] - O[idim]) / S[idim]) + 1
            Imax[idim] = ceil(Int64, (tmax[idim] - O[idim]) / S[idim]) + 1
        end

        # For raster indices falling into cell coordinate window,
        # check if raster points are in the cell
        # If so, obtain P1 interpolated raster data and
        # assign them to rasterflux
        for I[1] in Imin[1]:Imax[1] # 0 alloc
            if I[1] <= 0 || I[1] > length(rastercoord[1])
                continue
            end
            for I[2] in Imin[2]:Imax[2]
                if I[2] <= 0 || I[2] > length(rastercoord[2])
                    continue
                end
                for I[3] in Imin[3]:Imax[3]
                    if dim == 3 && (I[3] <= 0 || I[3] > length(rastercoord[3]))
                        continue
                    end

                    # Fill raster point to be tested
                    for idim in 1:dim
                        X[idim] = rastercoord[idim][I[idim]]
                    end

                    # Get barycentric coordinates
                    bary!(λ, invA, L2G, X)

                    # Check positivity of bc coordinates with some slack for
                    # round-off errors. Therefore a point may be found in two
                    # neighboring triangles. Constraining points to the raster ensures
                    # that only the last of them is taken.
                    if all(x -> x > -tol, λ)
                        # Interpolate vector value
                        if typeof(v) <: Function
                            v(V, λ, icell)
                        else
                            fill!(V, 0.0)
                            for inode in 1:(dim + 1)
                                for idim in 1:dim
                                    V[idim] += λ[inode] * v[idim, cn[inode, icell]]
                                end
                            end
                        end

                        for idim in 1:dim
                            rasterflux[idim, I[1], I[2], I[3]] = V[idim]
                        end
                    end
                end
            end
        end
    end
    return rastercoord, rasterflux
end

# Calculate isolevel values and function limits
function isolevels(ctx, funcs)
    return makeisolevels(
        funcs,
        ctx[:levels],
        ctx[:limits] == :auto ? (1, -1) : ctx[:limits],
        ctx[:colorbarticks] == :default ? nothing : ctx[:colorbarticks]
    )
end

"""
     function quiverdata(rastercoord,
                         rasterflux;
                         vscale = 1.0,
                         vnormalize = true,
                         vconstant = false)



Extract  nonzero fluxes for quiver plots from rastergrid obtained by [`vectorsample`](@ref).

Returns qc, qv -  `d x nquiver` matrices.

If vnormalize is true, the vector field is normalized to vscale*min(spacing), otherwise, it
is scaled by vscale
Result data are meant to  be ready for being passed to calls to `quiver` with various plotting backends.

"""
function quiverdata(
        rastercoord,
        rasterflux;
        vscale = 1.0,
        vnormalize = true,
        vconstant = false
    )
    dim = length(rastercoord)

    imax = length(rastercoord[1])
    jmax = length(rastercoord[2])
    spacing = (rastercoord[1][2] - rastercoord[1][1], rastercoord[2][2] - rastercoord[2][1])
    kmax = 1
    if dim > 2
        spacing = (spacing..., rastercoord[3][2] - rastercoord[3][1])
        kmax = length(rastercoord[3])
    end

    # memory for point, vector to be investigated
    X = zeros(Float32, dim)
    V = zeros(Float32, dim)

    qvcoord = Vector{SVector{dim, Float32}}(undef, 0)
    qvvalues = Vector{SVector{dim, Float32}}(undef, 0)

    I = ones(Int, 3)
    for I[1] in 1:imax # 0 allocs besides of push
        for I[2] in 1:jmax
            for I[3] in 1:kmax
                for idim in 1:dim
                    V[idim] = rasterflux[idim, I[1], I[2], I[3]]
                    X[idim] = rastercoord[idim][I[idim]]
                end
                if !iszero(V)
                    push!(qvcoord, X)
                    push!(qvvalues, V)
                end
            end
        end
    end

    # Reshape into matrices
    qc = reshape(reinterpret(Float32, qvcoord), (2, length(qvcoord)))
    qv = reshape(reinterpret(Float32, qvvalues), (2, length(qvvalues)))

    # Normalize vectors to raster point spacing
    if vconstant
        for j in 1:size(qv, 2)
            view(qv, :, j) ./= norm(view(qv, :, j))
        end
    elseif vnormalize
        @views vmax = maximum(norm, (qv[:, i] for i in 1:length(qvvalues)))
        vscale = vscale * min(spacing...) / vmax
    end

    # Scale vectors with user input
    qv .*= vscale
    return qc, qv
end

function bary!(λ, invA, L2G, x)
    mapderiv!(invA, L2G, nothing)
    fill!(λ, 0)
    for j in 1:length(x)
        dj = x[j] - L2G.b[j]
        for k in 1:length(x)
            λ[k] += invA[j, k] * dj
        end
    end
    ExtendableGrids.postprocess_xreftest!(λ, Triangle2D)
    return nothing
end

function cellcolors(grid, coloring)
    xr = grid[CellRegions]
    if coloring == :partitions
        xr = similar(xr)
        for icol in pcolors(grid)
            for ipart in pcolor_partitions(grid, icol)
                for icell in partition_cells(grid, ipart)
                    xr[icell] = ipart
                end
            end
        end
    elseif coloring == :pcolors
        xr = similar(xr)
        for icol in pcolors(grid)
            for ipart in pcolor_partitions(grid, icol)
                for icell in partition_cells(grid, ipart)
                    xr[icell] = icol
                end
            end
        end
    end
    return xr
end

function num_cellcolors(grid, coloring)
    if coloring == :partitions
        return num_partitions(grid)
    elseif coloring == :pcolors
        return num_pcolors(grid)
    else
        return num_cellregions(grid)
    end
end


const leglocs = Dict(
    :none => "",
    :best => "best",
    :lt => "upper left",
    :ct => "upper center",
    :rt => "upper right",
    :lc => "center left",
    :cc => "center",
    :rc => "center right",
    :lb => "lower left",
    :cb => "lower center",
    :rb => "lower right"
)
