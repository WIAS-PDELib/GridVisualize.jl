using Documenter, ExtendableGrids, Literate, GridVisualize, Pluto
using GridVisualize.FlippableLayout
import CairoMakie
CairoMakie.activate!(; type = "svg", visible = false)
using Test

plotting = joinpath(@__DIR__, "..", "examples", "plotting.jl")
include(plotting)

function rendernotebook(name)
    ENV["PLUTO_PROJECT"] = @__DIR__
    println("rendernotebook($(name))")
    input = joinpath(@__DIR__, "..", "examples", name * ".jl")
    output = joinpath(@__DIR__, "src", name * ".html")
    session = Pluto.ServerSession()
    notebook = Pluto.SessionActions.open(session, input; run_async = false)
    html_contents = Pluto.generate_html(notebook)
    return write(output, html_contents)
end

include("makeplots.jl")

example_md_dir = joinpath(@__DIR__, "src", "examples")

function mkdocs()
    Literate.markdown(plotting, example_md_dir; documenter = false, info = false)

    if true
        rendernotebook("plutovista")
        makeplots(example_md_dir; Plotter = CairoMakie, extension = "svg")
    end
    generated_examples = joinpath.("examples", filter(x -> endswith(x, ".md"), readdir(example_md_dir)))
    return makedocs(;
        sitename = "GridVisualize.jl",
        modules = [GridVisualize],
        doctest = false,
        clean = false,
        authors = "J. Fuhrmann",
        repo = "https://github.com/WIAS-PDELib/GridVisualize.jl",
        pages = [
            "Home" => "index.md",
            "Public API" => "api.md",
            "Private API" => "privapi.md",
            "Examples" => generated_examples,
        ]
    )
end

mkdocs()

deploydocs(; repo = "github.com/WIAS-PDELib/GridVisualize.jl.git", devbranch = "main")
