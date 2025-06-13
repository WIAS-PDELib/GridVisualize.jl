# # Plotting TriangulateIO
#
# GridVisualize comes with [`GridVisualize.plot_triangulateio`](@ref), a method to plot the input/output
# struct for Triangle mesh generator, provided by the [Triangulate.jl](https://github.com/JuliaGeometry/Triangulate.jl) wrapper.
# Supported are PyPlot and Makie backends.
#
#
using GridVisualize
using Triangulate
using Printf

function plotting_triangulateio(; Plotter = default_plotter(), maxarea = 0.05, resolution = (600, 300))
    triin = Triangulate.TriangulateIO()
    triin.pointlist = Matrix{Cdouble}([0.0 0.0; 1.0 0.0; 1.0 1.0; 0.6 0.6; 0.0 1.0]')
    triin.segmentlist = Matrix{Cint}([1 2; 2 3; 3 4; 4 5; 5 1]')
    triin.segmentmarkerlist = Vector{Int32}([1, 2, 3, 4, 5])
    area = @sprintf("%.15f", maxarea)
    (triout, vorout) = triangulate("pa$(area)DQ", triin)

    vis = GridVisualizer(;
        Plotter, layout = (1, 2), clear = true, resolution
    )

    plot_triangulateio!(vis[1, 1], triin, title = "input")
    plot_triangulateio!(
        vis[1, 2],
        triout;
        voronoi = vorout,
        circumcircles = true,
        title = "output"
    )
    return reveal(vis)
end;

# ![](plotting_triangulateio.png)


plotting_functions_png = [
    :plotting_triangulateio,
]

function generateplots(picdir; Plotter = nothing)
    filepaths = String[]
    if isdefined(Plotter, :Makie)
        size = (600, 300)
        Plotter.activate!(; type = "png", visible = false)
        for plotting_f in plotting_functions_png
            @eval begin
                path = joinpath($picdir, "$($plotting_f).png")
                p = $plotting_f(; Plotter = $Plotter)
                $Plotter.save(path, p)
                println("successfully generated plot for $($plotting_f)")
                push!($filepaths, path)
            end
        end
    end
    return filepaths
end
