using Test, ExtendableGrids, GridVisualize, Pkg, LinearAlgebra

import CairoMakie, PyPlot, PlutoVista

CairoMakie.activate!(; type = "svg", visible = false)

plotting = joinpath(@__DIR__, "..", "examples", "plotting.jl")
include(plotting)

for Plotter in [CairoMakie]
    @eval begin
        @testset "generateplots - $(nameof($Plotter))" begin
            mktempdir() do dir
                filepaths = generateplots(dir; Plotter = $Plotter)
                for path in filepaths
                    @test isfile(path)
                end
            end
        end
    end
end

# Some Plotters cannot perform the `makeplots` run, only try a `multiscene`
for Plotter in [PyPlot, PlutoVista]
    @eval begin
        @testset "plotting_multiscene - $(nameof($Plotter))" begin
            @test plotting_multiscene(Plotter = $Plotter) !== nothing
        end
    end
end

include("test_slice_plots.jl")


if isdefined(Docs, :undocumented_names) # >=1.11
    @testset "UndocumentedNames" begin
        @test isempty(Docs.undocumented_names(GridVisualize))
    end
end
