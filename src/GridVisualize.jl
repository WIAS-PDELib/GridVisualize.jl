"""
    GridVisualize

$(read(joinpath(@__DIR__, "..", "README.md"), String))
"""
module GridVisualize

using Printf
using LinearAlgebra

using DocStringExtensions
using OrderedCollections
using ElasticArrays
using StaticArrays
using Colors
using ColorSchemes
using GeometryBasics
using GridVisualizeTools
using ExtendableGrids

include("dispatch.jl")
include("common.jl")
include("slice_plots.jl")

export scalarplot, scalarplot!
export gridplot, gridplot!
export vectorplot, vectorplot!
export streamplot, streamplot!
export customplot, customplot!
export quiverdata, vectorsample
export plot_triangulateio, plot_triangulateio!
export save, reveal
export isplots, isvtkview, ispyplot, ismakie, isplutovista
export GridVisualizer, SubVisualizer
export plottertype, available_kwargs
export default_plotter!, default_plotter
export PyPlotType, MakieType, PlotsType, VTKViewType, PlutoVistaType, MeshCatType
export movie

end
