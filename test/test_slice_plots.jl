@testset "eval_expr! test" begin

    function test_slice_expr(result, expr)
        vec = zero(result)
        GridVisualize.eval_slice_expr!(vec, expr)
        @test vec == result
        return nothing
    end

    function test_fail_slice_plot(n, expr)
        try
            GridVisualize.eval_slice_expr!(zeros(n), expr)
        catch e
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

    test_slice_expr([1, 0, 0, -3], :x => 3)
    test_slice_expr([0, 1, 0, 5.5], :y => -5.5)
    test_slice_expr([0, 0, 1, -1 // 2], :z => 1 // 2)

    test_fail_slice_plot(3, :(y^2)) # wrong operation
    test_fail_slice_plot(4, :(z * z - 4)) # wrong operation
    test_fail_slice_plot(3, :(x * y)) # wrong operation
    test_fail_slice_plot(3, :(w + x - 5)) # wrong symbol
    test_fail_slice_plot(3, :(z + x - 5)) # :z not allowed in 2D

    test_fail_slice_plot(4, :w => 3) # wrong symbol
    test_fail_slice_plot(4, :x => "wtf") # wrong type
    test_fail_slice_plot(3, :z => 1) # :z not allowed in 2D
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


@testset "slice_plot!" begin

    let
        grid = simplexgrid(0.0:10.0, 0.0:10.0)
        f = map(x -> x[1] * x[2], grid)

        @test scalarplot(grid, f, Plotter = CairoMakie, slice = :(3x - y + 5)) !== nothing
        @test scalarplot(grid, f, Plotter = CairoMakie, slice = :y => 4) !== nothing
    end

    let
        grid = simplexgrid(0.0:10.0, 0.0:10.0, 0.0:10.0)
        f = map(x -> x[1] * x[2] + x[3], grid)

        @test scalarplot(grid, f, Plotter = CairoMakie, slice = :(3x - y - 5.5z + 5)) !== nothing
        @test scalarplot(grid, f, Plotter = CairoMakie, slice = :z => 4) !== nothing
    end
end
