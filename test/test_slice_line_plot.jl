@testset "eval_expr! test" begin

    function test_slice_expr(result, expr)
        vec = zero(result)
        GridVisualize.eval_slice_expr!(vec, expr)
        @test vec == result
        return nothing
    end

    function test_fail_slice_plot(expr)
        try
            GridVisualize.eval_slice_expr!(zeros(4), expr)
        catch
            @show "fail with $expr"
            @test true
        else
            @show expr
            @test false
        end
        return nothing
    end

    test_slice_expr([1, 0, 0, 0], :(x))
    test_slice_expr([0, 2, 0, 0], :(2y))
    test_slice_expr([0, 0, -3, 0], :(-3z))
    test_slice_expr([1, 2, 3, 4], :(3z + x + 2y + 4))
    test_slice_expr([-4, 3, -2, 1], :(1 - 2z + 3y - 4x))

    test_fail_slice_plot(:(y^2)) # wrong operation
    test_fail_slice_plot(:(z * z - 4)) # wrong operation
    test_fail_slice_plot(:(x * y)) # wrong operation
end


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
