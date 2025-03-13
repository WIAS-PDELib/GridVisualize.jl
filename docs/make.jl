using Documenter, ExtendableGrids, GridVisualize
using GridVisualize.FlippableLayout
import PlutoSliderServer # LoadError: Please import/use PlutoSliderServer.jl in order to use docplutonotebooks with `iframe=true`
using ExampleJuggler
import CairoMakie
ExampleJuggler.verbose!(true)
using Test

function mkdocs()
    cleanexamples()
    example_dir = joinpath(@__DIR__, "..", "examples")
    notebook_dir = joinpath(@__DIR__, "..", "examples")

    generated_examples = @docscripts(example_dir, ["Plotting Examples" => "plotting.jl"], Plotter = CairoMakie)
    notebook_examples = @docplutonotebooks(notebook_dir, ["plutovista.jl"], iframe = true, iframe_height = "2000px")

    makedocs(;
        sitename = "GridVisualize.jl",
        modules = [GridVisualize],
        doctest = false,
        clean = false,
        authors = "J. Fuhrmann",
        repo = "https://github.com/WIAS-PDELib/GridVisualize.jl",
        pages = [
            "Home" => "index.md",
            "Public API" => "api.md",
            "Examples" => [
                "Plotting examples" => generated_examples,
                "Pluto notebooks" => notebook_examples,
            ],
            "Private API" => "privapi.md",
            "Contributing" => "contributing.md",
        ]
    )
    return nothing
end

mkdocs()

deploydocs(; repo = "github.com/WIAS-PDELib/GridVisualize.jl.git", devbranch = "main")
