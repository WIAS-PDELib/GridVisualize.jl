import VTKView
import Plots
import PyPlot
import GLMakie

include("examples/plotting.jl")

function makefigs()
    Plotter = PyPlot
    scene = plotting_multiscene(; Plotter = Plotter)
    save("multiscene_pyplot.png", scene; Plotter = Plotter)

    Plotter = Plots
    scene = plotting_multiscene(; Plotter = Plotter)
    save("multiscene_plots.png", scene; Plotter = Plotter)

    Plotter = GLMakie
    scene = plotting_multiscene(; Plotter = Plotter)
    save("multiscene_glmakie.png", scene; Plotter = Plotter)

    Plotter = VTKView
    scene = plotting_multiscene(; Plotter = Plotter)
    return save("multiscene_vtkview.png", scene; Plotter = Plotter)
end
