using Documenter, GridVisualize
import PlutoSliderServer # LoadError: Please import/use PlutoSliderServer.jl in order to use docplutonotebooks with `iframe=true`
using ExampleJuggler
import CairoMakie, PlutoVista, MeshCat, VTKView # for docstrings and graphics generated by documentet
import Plots, PyPlot # for docstrings only
using ExtendableGrids
ExampleJuggler.verbose!(true)
using Test

GridVisualizeMakieExt = Base.get_extension(GridVisualize, :GridVisualizeMakieExt)
GridVisualizePlotsExt = Base.get_extension(GridVisualize, :GridVisualizePlotsExt)
GridVisualizePlutoVistaExt = Base.get_extension(GridVisualize, :GridVisualizePlutoVistaExt)
GridVisualizePyPlotExt = Base.get_extension(GridVisualize, :GridVisualizePyPlotExt)
GridVisualizeMeshCatExt = Base.get_extension(GridVisualize, :GridVisualizeMeshCatExt)
GridVisualizeVTKViewExt = Base.get_extension(GridVisualize, :GridVisualizeVTKViewExt)


function mkdocs()
    cleanexamples()
    example_dir = joinpath(@__DIR__, "..", "examples")
    notebook_dir = joinpath(@__DIR__, "..", "examples")
    generated_examples = @docscripts(
        example_dir, [
            "Plotting Examples" => "plotting.jl",
        ], Plotter = CairoMakie
    )
    notebook_examples = @docplutonotebooks(notebook_dir, ["plutovista.jl"], iframe = true, iframe_height = "2000px")
    makedocs(;
        sitename = "GridVisualize.jl",
        modules = [
            GridVisualize,
            GridVisualizeMakieExt,
            GridVisualizePlutoVistaExt,
            GridVisualizePyPlotExt,
            GridVisualizePlotsExt,
            GridVisualizeVTKViewExt,
            GridVisualizeMeshCatExt,
        ],
        doctest = false,
        clean = false,
        authors = "J. Fuhrmann and contributors",
        repo = "https://github.com/WIAS-PDELib/GridVisualize.jl",
        pages = [
            "Home" => "index.md",
            "Changelog" => "changes.md",
            "Public API" => "api.md",
            "Examples" => generated_examples,
            "Pluto notebooks" => notebook_examples,
            "Private API" => "privapi.md",
            "Contributing" => "contributing.md",
        ]
    )
    return nothing

end

mkdocs()

deploydocs(; repo = "github.com/WIAS-PDELib/GridVisualize.jl.git", devbranch = "main")
