# Changelog

## [unreleased]

### Added

- `gridplot` with `GLMakie` and `PyPlot` now annotates the color bars with `boundary regions` and `cell regions`

- unit test for undocumented names

### Removed

- dangling exported function `backend!`

## [1.9.0] - 2024-11-25

### Added

- new key word argument `show_colorbar` (default=true) to toggle the color bar next to grid plots

### Fixed

- removed dependency on type piracy of Colors.RGB

## [1.8.2] - 2024-11-08

### Fixed

- Fix GLMakie region colorbar in 2D/3D Plots (#23)

## [1.8.1] - 2024-11-02
- update links after move to WIAS-PDELib org

## [1.8.0] - 2024-09-29

- Allow plotting using grids without boundary info

- Add dispatches for plotting with 'coord, cellnodes'


## [1.7.0] - 2024-06-18

- Add :cellcoloring keyword for handling partitioned grids

- Upgrade project.toml, dependencies

- Fix cellcolor numbering

- Remove nightly from ci due to pluto problem

- Merge pull request #27 from j-fu/handle-partitioning

Handle partitioning

## [1.6.0] - 2024-05-24

* require Julia >= 1.9

* allow for Makie 0.21;

See also https://blog.makie.org/blogposts/v0.21/



## [1.5.0] - 2023-12-09

- Gridscale etc (#20)

Some updates, fixes:

* gridscale for plutovista, pyplot, makie,plots

* Export vectorsample and quiverdata

* Fix streamplot handling for Makie

* spacing -> rasterpoints for quiver, streamplot

* Ensure that colorbarticks are always shown and contain the function limits

* Add customplot



## [1.4.0] - 2023-12-05

- Prevent PyPlot from normalizing quiver vectors

- Pin CairoMakie version due to https://github.com/MakieOrg/Makie.jl/issues/3440

- Add warnings for functionality not implemented in Plots

- Bump CairoMakie dependency

- Add streamplot for Makie

- Remove another .px_area

- Update multiscene plot for makie



## [1.3.0] - 2023-11-15

- Update Makie, GridVisualizeTools versions

## [1.2.0] - 2023-11-11

- Fix notebook html generation

- Bgcolor->backgroundcolor in makie.jl

- Add compat for stdlibs, bump minor version


## [1.1.7] - 2023-10-12

- Sets default value false for kwarg vconstant introduced in 1.15 (#18)

Co-authored-by: Christian Merdon <merdon@wias-berlin.de>
- Fix Documenter v1 issues


## [1.1.5] - 2023-09-11

- Allow for PlutVista 1.0

- Improved PyPlot backend: (#17)

* better fontsize recognition

* correct fig sizes

* tight_layout() also for SubVisualizer reveal

* fixed rare clipping of last colorlevel in scalarplot

* coordinate limits (xlimits etc.) are used in vectorsample (such that scaling only is applied to vectors in the clipped area)

* new vector scaling method vconstant that scales all arrows equally

* repaired streamplot (U and V arguments needed transposition)

* added density argument to streamplot



## [1.1.0] - 2023-06-02

- Add weakdeps + compat for Makie & co to Project.toml

- Enable levelalpha/planealpha for makie
fix colorbarticks

- Try to fix isoline rendering in makie 2d



## [1.0.0] - 2023-02-05

- Subgridplots (#16)

Handle plots of discontinuous functions in Makie,Pyplot, PlutoVista
