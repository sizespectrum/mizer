# Plot the total yield of each species by gear through time

After running a projection, the total yield of each species by fishing
gear can be plotted against time.

## Usage

``` r
plotYieldGear(
  object,
  species = NULL,
  gears = NULL,
  total = FALSE,
  log_x = FALSE,
  log_y = TRUE,
  log = NULL,
  ylim = c(NA, NA),
  tlim = c(NA, NA),
  highlight = NULL,
  return_data = FALSE,
  ...
)
```

## Arguments

- object:

  An object of class
  [MizerSim](https://sizespectrum.org/mizer/reference/MizerSim-class.md)

- species:

  The species to be selected. Optional. By default all target species
  are selected. A vector of species names, or a numeric vector with the
  species indices, or a logical vector indicating for each species
  whether it is to be selected (TRUE) or not.

- gears:

  A vector of gear names to be included in the plot. Default is all
  gears.

- total:

  A boolean value that determines whether the total yield from all
  species is plotted as well. Default is FALSE.

- log_x:

  If `TRUE`, use a log10 x-axis. Default is `FALSE`.

- log_y:

  If `TRUE`, use a log10 y-axis. Default is `TRUE`.

- log:

  Character string specifying which axes should use log10 scales, in the
  same form as the base
  [`plot()`](https://sizespectrum.org/mizer/reference/plot.md) argument.
  For example, `"x"`, `"y"`, `"xy"` or `""`. If supplied, this overrides
  `log_x` and `log_y`. For backward compatibility, `TRUE` and `FALSE`
  are interpreted as setting only `log_y`.

- ylim:

  A numeric vector of length two providing lower and upper limits for
  the y axis. Use `NA` to refer to the existing minimum or maximum.

- tlim:

  A numeric vector of length two providing lower and upper limits for
  the time axis, e.g. `c(1980, 2000)`. Use `NA` to apply no limit at
  that end. Default is `c(NA, NA)`.

- highlight:

  Name or vector of names of the species to be highlighted.

- return_data:

  A boolean value that determines whether the formatted data used for
  the plot is returned instead of the plot itself. Default is FALSE.

- ...:

  Arguments passed to
  [`getYieldGear()`](https://sizespectrum.org/mizer/reference/getYieldGear.md).

## Value

A ggplot2 object, unless `return_data = TRUE`, in which case a data
frame with the four variables 'Year', 'Yield', 'Species' and 'Gear' is
returned.

## Details

This plot is pretty easy to do by hand. It just gets the biomass using
the
[`getYieldGear()`](https://sizespectrum.org/mizer/reference/getYieldGear.md)
method and plots using the ggplot2 package. You can then fiddle about
with colours and linetypes etc. Just look at the source code for
details.

## See also

[plotting_functions](https://sizespectrum.org/mizer/reference/plotting_functions.md),
[`getYieldGear()`](https://sizespectrum.org/mizer/reference/getYieldGear.md)

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
[`plotSpectra2()`](https://sizespectrum.org/mizer/reference/plotSpectra2.md),
[`plotSpectraRelative()`](https://sizespectrum.org/mizer/reference/plotSpectraRelative.md),
[`plotYield()`](https://sizespectrum.org/mizer/reference/plotYield.md),
[`plotting_functions`](https://sizespectrum.org/mizer/reference/plotting_functions.md)

## Examples

``` r
# \donttest{
params <-  NS_params
sim <- project(params, effort=1, t_max=20, t_save = 0.2, progress_bar = FALSE)
plotYieldGear(sim)

plotYieldGear(sim, species = c("Cod", "Herring"), total = TRUE)


# Returning the data frame
fr <- plotYieldGear(sim, return_data = TRUE)
str(fr)
#> 'data.frame':    1212 obs. of  4 variables:
#>  $ Year   : num  0 0.2 0.4 0.6 0.8 1 1.2 1.4 1.6 1.8 ...
#>  $ Yield  : num  2.21e+11 1.86e+11 1.60e+11 1.42e+11 1.30e+11 ...
#>  $ Species: Factor w/ 12 levels "Sprat","Sandeel",..: 1 1 1 1 1 1 1 1 1 1 ...
#>  $ Gear   : Factor w/ 4 levels "Industrial","Pelagic",..: 1 1 1 1 1 1 1 1 1 1 ...
# }
```
