using Test, ExtendableGrids, GridVisualize, Pkg
import CairoMakie
CairoMakie.activate!(; type = "svg", visible = false)

plotting = joinpath(@__DIR__, "..", "examples", "plotting.jl")
include(plotting)
include("../docs/makeplots.jl")
@testset "makeplots - CairoMakie" begin
    makeplots(mktempdir(); Plotter = CairoMakie, extension = ".svg")
end


if isdefined(Docs, :undocumented_names) # >=1.11
    @testset "UndocumentedNames" begin
        @test isempty(Docs.undocumented_names(GridVisualize))
    end
end
