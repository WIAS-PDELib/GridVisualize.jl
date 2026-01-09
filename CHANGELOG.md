# Changelog

## [1.17.0]
- Add `:symlog` axis scaling for symmetric logarithmic axis scaling, see, e.g. [matplotlib symlog](https://matplotlib.org/stable/gallery/scales/symlog_demo.html)

## [1.16.3] - 2026-01-07
- Add `title`, `xlabel`, and `ylabel` to 1D and 2D `gridplot!` with `GLMakie`
- Add `title`, `xlabel`, and `ylabel` to 1D and 2D `gridplot!` with `Py[thon]Plot`

## [1.16.2] - 2025-12-04
- `leglocs[:cc]`: correct "center center" to "center"

## [1.16.1] - 2025-12-01
- `Plots.jl` and `PyPlot.jl`: fixed missing plots in VSCode by using `display` in `reveal`

## [1.16.0] - 2025-11-18
- Add PythonPlot to the backends

## [1.15.3] - 2025-06-24
- Fixed bug with GLMakie based animated plots from the REPL
- Added examples for animated plots from the REPL

## [1.15.2] - 2025-06-24
- fixed a corner case in `compute_2d_rotation_matrix` leading to wrong coordinate transformation in 1D slice plots

## [1.15.1] - 2025-06-22
- Allow for Makie 0.24 and Triangulate 3.0

## [1.15.0] - 2025-06-17
- Removed automatic overrides of user-supplied keywords `clear` , `show`, and `reveal` in `scalarplot!`, `vectorplot!`, `streamplot!`, and `customplot!`.
- Included `CHANGELOG.md` in docs.

## [1.14.0] - 2025-06-13
- Move plotting code to extensions
- Add `plot_triangulateio` method moved from [Triangulate.jl](https://github.com/JuliaGeometry/Triangulate.jl)

## [1.13.0] - 2025-06-12
- Bump makie dependency

## [1.12.0] - 2025-05-12
- Manage documentation examples with ExampleJuggler

## [1.11.0] - 2025-03-04

### Added

 - New key word argument `slice` for `scalarplot` to extract a `dim-1` scalar plot along a slice definition.
   The slice can be provided as an expression `slice = :(x + 2y - 4)` constrained by zero, or as a pair `slice = :x => 3`
   to fix a certain value of one axis.

## [1.10.0] - 2025-01-20

### Added

- `gridplot` with `GLMakie` and `PyPlot` now annotates the color bars with `boundary regions` and `cell regions`

- unit test for undocumented names

- Unit test for `PyPlot` and `PlutoVista`

- Allow to use breaking releases GeometryBasics v0.5 and Makie v0.22

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
