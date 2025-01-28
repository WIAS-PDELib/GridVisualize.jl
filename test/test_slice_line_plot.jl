@testset "z-rotation matrix" begin

    # t is already the z-vector ⇒ R = id
    t = [0, 0, 1]
    @test GridVisualize.compute_3d_z_rotation_matrix(t) ≈ I

    # random test
    t = [0.3543579344795591, 0.21250427156069385, 0.4465606445966942]
    R = GridVisualize.compute_3d_z_rotation_matrix(t)

    ξ = R[:, 1]
    η = R[:, 2]

    t ./= norm(t)
    @test R[:, 3] ≈ t

    # check that new coordinate system is orthogonal
    @test abs(ξ ⋅ η) < 1.0e-15
    @test abs(η ⋅ t) < 1.0e-15
    @test abs(ξ ⋅ t) < 1.0e-15

    # corner case I: only one axis rotation necessary
    t = [1, 0, 0]
    @test GridVisualize.compute_3d_z_rotation_matrix(t) ≈ [0 0 1; 0 1 0; -1 0 0]

    # corner case II: only one axis rotation necessary
    t = [0, 1, 0]
    @test GridVisualize.compute_3d_z_rotation_matrix(t) ≈ [1 0 0; 0 0 1; 0 -1 0]

end


@testset "2d-rotation matrix" begin

    # t is already the y-vector ⇒ R = id
    t = [0, 1]
    @test GridVisualize.compute_2d_rotation_matrix(t) ≈ I

    # random test
    t = [0.3543579344795591, 0.21250427156069385]
    R = GridVisualize.compute_2d_rotation_matrix(t)

    # test rotation property
    @test R't ≈ [0, 1]
end


@testset "lineplot" begin

    grid = simplexgrid(0.0:10.0, 0.0:10.0)

    f = map(x -> x[1] * x[2], grid)

    # x = 5 line
    line = [1, 0, -5]

    @test lineplot(grid, f, line, Plotter = CairoMakie) !== nothing
end


@testset "sliceplot" begin

    grid = simplexgrid(0.0:10.0, 0.0:10.0, 0.0:10.0)

    f = map(x -> x[1] * x[2] + x[3], grid)

    # x = 5 line
    plane = [1, 0, 0, -5]

    @test sliceplot(grid, f, plane, Plotter = CairoMakie) !== nothing
end
