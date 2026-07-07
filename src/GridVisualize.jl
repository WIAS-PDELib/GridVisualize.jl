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
include("pycommon.jl")
include("slice_plots.jl")
include("deprecated.jl")

export scalarplot, scalarplot!
export gridplot, gridplot!
export vectorplot, vectorplot!
export streamplot, streamplot!
export customplot, customplot!
export quiverdata, vectorsample
export plot_triangulateio, plot_triangulateio!
export save, reveal
export GridVisualizer, SubVisualizer
export plottertype, available_kwargs
export default_plotter!, default_plotter
export PlotterType, PyPlotType, PythonPlotType, UnionPythonPlotterType, UnionMakieType, PlotsType, VTKViewType, PlutoVistaType, MeshCatType, UnicodePlotsType
export isvtkview, ispyplot, ispythonplot, isplots, ismakie, ismeshcat, isplutovista, isunicodeplots
export movie

end
