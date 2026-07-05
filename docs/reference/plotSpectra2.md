# Compare abundance and biomass spectra from two objects

`plotSpectra2()` compares the abundance spectra from two `MizerParams`
or `MizerSim` objects in a single plot. Colours identify species or
groups and linetype identifies the object.

## Usage

``` r
plotSpectra2(
  object1,
  object2,
  name1 = "First",
  name2 = "Second",
  species = NULL,
  wlim = c(NA, NA),
  llim = c(NA, NA),
  ylim = c(NA, NA),
  power = 1,
  biomass = TRUE,
  total = FALSE,
  resource = TRUE,
  background = TRUE,
  highlight = NULL,
  log_x = TRUE,
  log_y = TRUE,
  log = NULL,
  size_axis = c("w", "l"),
  ...
)
```

## Arguments

- object1:

  First `MizerParams` or `MizerSim` object.

- object2:

  Second `MizerParams` or `MizerSim` object.

- name1, name2:

  Labels for the two objects, used in the linetype legend.

- species:

  The species to be selected. Optional. By default all target species
  are selected. A vector of species names, or a numeric vector with the
  species indices, or a logical vector indicating for each species
  whether it is to be selected (TRUE) or not.

- wlim:

  A numeric vector of length two providing lower and upper limits for
  the w axis. Use NA for the default: the lower default is
  `min(params@w) / 100` when `resource = TRUE` (to show some resource
  below the fish grid) or `min(params@w)` when `resource = FALSE`; the
  upper default is `max(params@w_full)`. Data is filtered to this range
  and the axis limits are set accordingly.

- llim:

  A numeric vector of length two providing lower and upper limits for
  the length axis when `size_axis = "l"`. Use `NA` to auto-scale to the
  data range. Data is filtered to this range and the axis limits are set
  accordingly.

- ylim:

  A numeric vector of length two providing lower and upper limits for
  the y axis. Use NA to auto-scale to the data range. Values below 1e-20
  are always filtered out from the data regardless of `ylim[1]`. Data
  above `ylim[2]` is filtered and the upper axis limit is set
  accordingly.

- power:

  The abundance is plotted as the number density times the weight raised
  to `power`. The default `power = 1` gives the biomass density, whereas
  `power = 2` gives the biomass density with respect to logarithmic size
  bins.

- biomass:

  **\[deprecated\]** Only used if `power` argument is missing. Then
  `biomass = TRUE` is equivalent to `power=1` and `biomass = FALSE` is
  equivalent to `power=0`

- total:

  A boolean value that determines whether the total over all species in
  the system is plotted as well. Note that even if the plot only shows a
  selection of species, the total is including all species. Default is
  FALSE.

- resource:

  A boolean value that determines whether resource is included. Default
  is TRUE.

- background:

  A boolean value that determines whether background species are
  included. Ignored if the model does not contain background species.
  Default is TRUE.

- highlight:

  Name or vector of names of the species to be highlighted by being
  plotted with thicker lines.

- log_x:

  If `TRUE` (default), use a log10 x-axis.

- log_y:

  If `TRUE` (default), use a log10 y-axis.

- log:

  Character string specifying which axes should use log10 scales, in the
  same form as the base
  [`plot()`](https://sizespectrum.org/mizer/reference/plot.md) argument.
  For example, `"x"`, `"y"`, `"xy"` or `""`. If supplied, this overrides
  `log_x` and `log_y`.

- size_axis:

  Whether to plot size as weight (`"w"`, default) or length (`"l"`),
  using the allometric weight-length relationship.

- ...:

  Additional arguments passed to
  [`plotSpectra()`](https://sizespectrum.org/mizer/reference/plotSpectra.md)
  for preparing the spectra data, for example `time_range` or
  `geometric_mean` for `MizerSim` objects.

## Value

A ggplot2 object. `plotlySpectra2()` returns a plotly object.

## Details

`plotlySpectra2()` is the interactive plotly version.

## See also

[plotting_functions](https://sizespectrum.org/mizer/reference/plotting_functions.md),
[`plotSpectra()`](https://sizespectrum.org/mizer/reference/plotSpectra.md),
[`plotSpectraRelative()`](https://sizespectrum.org/mizer/reference/plotSpectraRelative.md)

Other plotting functions:
[`addPlot()`](https://sizespectrum.org/mizer/reference/addPlot.md),
[`animate.ArrayTimeBySpeciesBySize()`](https://sizespectrum.org/mizer/reference/animate.md),
[`plot`](https://sizespectrum.org/mizer/reference/plot.md),
[`plot2()`](https://sizespectrum.org/mizer/reference/plot2.md),
[`plotBiomass()`](https://sizespectrum.org/mizer/reference/plotBiomass.md),
[`plotCDF()`](https://sizespectrum.org/mizer/reference/plotCDF.md),
[`plotCDF2()`](https://sizespectrum.org/mizer/reference/plotCDF2.md),
[`plotDiet()`](https://sizespectrum.org/mizer/reference/plotDiet.md),
[`plotFMort()`](https://sizespectrum.org/mizer/reference/plotFMort.md),
[`plotFeedingLevel()`](https://sizespectrum.org/mizer/reference/plotFeedingLevel.md),
[`plotGrowthCurves()`](https://sizespectrum.org/mizer/reference/plotGrowthCurves.md),
[`plotMizerParams`](https://sizespectrum.org/mizer/reference/plotMizerParams.md),
[`plotMizerSim`](https://sizespectrum.org/mizer/reference/plotMizerSim.md),
[`plotPredMort()`](https://sizespectrum.org/mizer/reference/plotPredMort.md),
[`plotRelative()`](https://sizespectrum.org/mizer/reference/plotRelative.md),
[`plotSpectra()`](https://sizespectrum.org/mizer/reference/plotSpectra.md),
[`plotSpectraRelative()`](https://sizespectrum.org/mizer/reference/plotSpectraRelative.md),
[`plotYield()`](https://sizespectrum.org/mizer/reference/plotYield.md),
[`plotYieldGear()`](https://sizespectrum.org/mizer/reference/plotYieldGear.md),
[`plotting_functions`](https://sizespectrum.org/mizer/reference/plotting_functions.md)

## Examples

``` r
# \donttest{
sim1 <- project(NS_params, t_max = 10, progress_bar = FALSE)
sim2 <- project(NS_params, effort = 0.5, t_max = 10, progress_bar = FALSE)
plotSpectra2(sim1, sim2, "Original", "Effort = 0.5")

# }
```
