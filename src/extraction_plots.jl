"""
    $(SIGNATURES)

    Compute and return a rotation matrix 𝐑³ˣ³
    which rotates the third unit vector z = [0, 0, 1]ᵀ
    in the direction of a given target_vector t ∈ 𝐑³

    target_vector: vector with 3 components (not necessarily of unit length)
"""
function compute_3d_z_rotation_matrix(target_vector)

    # rotation by x-axis
    R_x(α) = @SArray [ 1 0 0 ; 0 cos(α) -sin(α); 0 sin(α) cos(α) ]

    # rotation by y-axis
    R_y(α) = @SArray [ cos(α)  0 sin(α); 0 1 0; -sin(α) 0 cos(α) ]

    # normalize t
    t = target_vector ./= norm(target_vector)

    # solve non-linear system
    α = asin(-t[2])

    if abs(α) ≈ π / 2 # cos(α) = 0: z → ±y, no y-rotation
        β = 0
    else
        β = acos(t[3] / cos(α))
    end

    return R_y(β) * R_x(α)
end


"""
    $(SIGNATURES)

    Compute and return a rotation matrix 𝐑²ˣ²
    which rotates the second unit vector y = [0, 1]ᵀ
    in the direction of a given target_vector t ∈ 𝐑²

    target_vector: vector with 2 components (not necessarily of unit length)
"""
function compute_2d_rotation_matrix(target_vector)

    # rotation matrix
    R(α) = @SArray [ cos(α) sin(α); -sin(α) cos(α) ]

    # normalize t
    t = target_vector ./= norm(target_vector)

    # solve non-linear system for α
    α = asin(t[1])

    return R(α)
end

"""
    $(SIGNATURES)

    Extract a 2D slice plot of a 3D plot

    The intersection of the given 3D grid with a given plane is computed and rotated into a 2D
    coordinate system.

    For a simple plane (x = const / y = const / z = const) the original coordinates of the free
    axes are preserved. The plotted axes order is in this case
        x = const: ( y, z )
        y = const: ( x, z )
        z = const: ( x, y )

    Else, for a generic plane, the new coordinate system has non-negative values and start at [0,0].

    vis:     GridVisualizer
    grid:    3D ExtendableGrid
    values:  value vector corresponding to the grid nodes
    plane:   Vector [a,b,c,d], s.t., ax + by + cz + d = 0 defines the plane that slices the 3D grid
    xlabel:  new first transformed coordinate
    ylabel:  new second transformed coordinate
"""
function sliceplot!(vis, grid, values, plane; xlabel = "ξ", ylabel = "η", kwargs...)

    @assert length(plane) == 4 "a plane equation requires exactly 4 parameters"

    # get new data from marching_tetrahedra
    new_coords, new_triangles, new_values = GridVisualize.marching_tetrahedra(grid, values, [plane], [])

    # construct new 2D grid
    grid_2d = ExtendableGrid{Float64, Int32}()

    a::Float64, b::Float64, c::Float64, _ = plane

    # transformation matrix
    if (a, b, c) == (1.0, 0.0, 0.0)
        rotation_matrix = @SArray [
            0.0 0.0 1.0
            1.0 0.0 0.0
            0.0 1.0 0.0
        ]
    elseif (a, b, c) == (0.0, 1.0, 0.0)
        rotation_matrix = @SArray [
            1.0 0.0 0.0
            0.0 0.0 1.0
            0.0 1.0 0.0
        ]
    elseif (a, b, c) == (0.0, 0.0, 1.0)
        rotation_matrix = @SArray [
            1.0 0.0 0.0
            0.0 1.0 0.0
            0.0 0.0 1.0
        ]
    else
        rotation_matrix = compute_3d_z_rotation_matrix([a, b, c])
    end

    grid_2d[Coordinates] = Matrix{Float64}(undef, 2, length(new_coords))
    for (ip, p) in enumerate(new_coords)
        # to obtain the projected coordinates, we can simply use the transpose of the rotation matrix
        @views grid_2d[Coordinates][:, ip] .= (rotation_matrix'p)[1:2]
    end

    if (a, b, c) ∉ [(1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)]
        # adjust the coordinates s.t. the minimal coordinate is zero
        @views grid_2d[Coordinates][1, :] .-= minimum(grid_2d[Coordinates][1, :])
        @views grid_2d[Coordinates][2, :] .-= minimum(grid_2d[Coordinates][2, :])
    end

    grid_2d[CellNodes] = Matrix{Int32}(undef, 3, length(new_triangles))
    for (it, t) in enumerate(new_triangles)
        @views grid_2d[CellNodes][:, it] .= t
    end

    scalarplot!(vis, grid_2d, new_values; xlabel, ylabel, kwargs...)

    return vis
end

"""
    ($SIGNATURES)

    Handy in-place variant of [`sliceplot!`](@ref)
"""
function sliceplot(grid, values, line; Plotter = default_plotter(), kwargs...)
    return sliceplot!(GridVisualizer(; Plotter = Plotter, show = true, kwargs...), grid, values, line; kwargs...)
end


"""
    $(SIGNATURES)

    Extract a 1D line plot from a 2D plot

    The intersection of the given 2D grid with a given line is computed and rotated onto the
    first coordinate axis.

    For a simple line (x = const / y = const) the original coordinates of the other axis are
    preserved.

    Else, for a generic line, the new axis has non-negative values and starts at [0,0].

    vis:     GridVisualizer
    grid:    2D ExtendableGrid
    values:  value vector corresponding to the grid nodes
    line:    Vector [a,b,c], s.t., ax + by + d = 0 defines the line that slices the 2D grid
    xlabel:  new coordinate of the resulting line
    ylabel:  label for the data
"""
function lineplot!(vis, grid, values, line; xlabel = "line", ylabel = "value", kwargs...)

    @assert length(line) == 3 "a line equation requires exactly 3 parameters"

    # get new data from marching_triangles
    new_coords, new_adjecencies, new_values = GridVisualize.marching_triangles(grid, values, [line], [])

    # construct new 1D grid
    grid_1d = ExtendableGrid{Float64, Int32}()

    a::Float64, b::Float64, _ = line

    if (a, b) == (1.0, 0.0)
        rotation_matrix = @SArray [
            0.0 1.0
            1.0 0.0
        ]
    elseif (a, b) == (0.0, 1.0)
        rotation_matrix = @SArray [
            1.0 0.0
            0.0 1.0
        ]
    else
        rotation_matrix = compute_2d_rotation_matrix([a, b])
    end

    # rotated coordinates (drop second component)
    rot_coords = [ (rotation_matrix'p)[1] for p in new_coords ]

    grid_1d[Coordinates] = Matrix{Float64}(undef, 1, length(rot_coords))
    for ip in eachindex(rot_coords)
        # to obtain the projected coordinates, we can simply use the transpose of the rotation matrix
        grid_1d[Coordinates][1, ip] = rot_coords[ip]
    end

    if (a, b) ∉ [(1.0, 0.0), (0.0, 1.0)]
        # adjust the coordinates s.t. the minimal coordinate is always zero
        grid_1d[Coordinates] .-= minimum(grid_1d[Coordinates])
    end

    grid_1d[CellNodes] = Matrix{Int32}(undef, 2, length(new_adjecencies))
    for (it, t) in enumerate(new_adjecencies)
        @views grid_1d[CellNodes][:, it] .= t
    end

    # sort the coordinates: this is currently required by the Makie backend
    p = @views sortperm(grid_1d[Coordinates][:])
    grid_1d[Coordinates] = @views grid_1d[Coordinates][:, p]
    new_values = new_values[p]

    scalarplot!(vis, grid_1d, new_values; xlabel, ylabel, kwargs...)

    return vis

end

"""
    ($SIGNATURES)

    Handy in-place variant of [`lineplot!`](@ref)
"""
function lineplot(grid, values, line; Plotter = default_plotter(), kwargs...)
    return lineplot!(GridVisualizer(; Plotter = Plotter, show = true, kwargs...), grid, values, line; kwargs...)
end
