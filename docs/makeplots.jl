function makeplots(picdir; Plotter = GLMakie, extension = "png")
    p = plotting_multiscene(; Plotter = Plotter)

    fname = joinpath(picdir, "plotting_multiscene." * extension)
    save(fname, p; Plotter = Plotter)
    @test isfile(fname)
    println("multiscene")

    p = plotting_func1d(; Plotter = Plotter)
    fname = joinpath(picdir, "plotting_func1d." * extension)
    save(fname, p; Plotter = Plotter)
    @test isfile(fname)
    println("func1d")

    p = plotting_func2d(; Plotter = Plotter)
    fname = joinpath(picdir, "plotting_func2d." * extension)
    save(fname, p; Plotter = Plotter)
    @test isfile(fname)
    println("func2d")

    p = plotting_func3d(; Plotter = Plotter)
    fname = joinpath(picdir, "plotting_func3d." * extension)
    save(fname, p; Plotter = Plotter)
    @test isfile(fname)
    println("func3d")

    p = plotting_slice3d(; Plotter = Plotter)
    fname = joinpath(picdir, "plotting_slice3d." * extension)
    save(fname, p; Plotter = Plotter)
    @test isfile(fname)
    println("slice3d")

    p = plotting_line2d(; Plotter = Plotter)
    fname = joinpath(picdir, "plotting_line2d." * extension)
    save(fname, p; Plotter = Plotter)
    @test isfile(fname)
    println("line2d")

    fname = joinpath(picdir, "plotting_jfunc1d." * "gif")
    p = plotting_jfunc1d(; Plotter = Plotter, filename = fname)
    @test isfile(fname)
    println("jfunc1d")

    p = plotting_jfunc2d(; Plotter = Plotter)
    fname = joinpath(picdir, "plotting_jfunc2d." * extension)
    save(fname, p; Plotter = Plotter)
    @test isfile(fname)
    println("jfunc2d")

    p = plotting_jfunc3d(; Plotter = Plotter)
    fname = joinpath(picdir, "plotting_jfunc3d." * extension)
    save(fname, p; Plotter = Plotter)
    @test isfile(fname)
    println("jfunc3d")

    p = plotting_vec2d(; Plotter = Plotter)
    fname = joinpath(picdir, "plotting_vec2d." * extension)
    save(fname, p; Plotter = Plotter)
    @test isfile(fname)

    println("vec2d")

    p = plotting_stream2d(; Plotter = Plotter)
    fname = joinpath(picdir, "plotting_stream2d." * extension)
    save(fname, p; Plotter = Plotter)
    @test isfile(fname)
    println("stream2d")

    p = plotting_grid1d(; Plotter = Plotter)
    fname = joinpath(picdir, "plotting_grid1d." * extension)
    save(fname, p; Plotter = Plotter)
    @test isfile(fname)
    println("grid1d")

    p = plotting_grid2d(; Plotter = Plotter)
    fname = joinpath(picdir, "plotting_grid2d." * extension)
    save(fname, p; Plotter = Plotter)
    @test isfile(fname)
    println("grid2d")

    p = plotting_grid3d(; Plotter = Plotter)
    fname = joinpath(picdir, "plotting_grid3d." * extension)
    save(fname, p; Plotter = Plotter)
    @test isfile(fname)
    println("grid3d")

    fname = joinpath(picdir, "plotting_movie." * "gif")
    p = plotting_movie(; filename = fname, Plotter = Plotter)
    @test isfile(fname)
    println("movie")

    p = plotting_custom(; Plotter = Plotter)
    fname = joinpath(picdir, "plotting_custom." * extension)
    save(fname, p; Plotter = Plotter)
    @test isfile(fname)
    println("custom")

    return true
end
