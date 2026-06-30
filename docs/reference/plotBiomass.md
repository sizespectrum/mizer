# Plot the biomass of species through time

After running a projection, the biomass of each species can be plotted
against time. The biomass is calculated within user defined size limits
(see
[`getBiomass()`](https://sizespectrum.org/mizer/reference/getBiomass.md)).

## Usage

``` r
plotBiomass(
  object,
  species = NULL,
  tlim = c(NA, NA),
  y_ticks = 6,
  ylim = c(NA, NA),
  total = FALSE,
  background = TRUE,
  highlight = NULL,
  log = NULL,
  return_data = FALSE,
  log_x = FALSE,
  log_y = TRUE,
  use_cutoff = FALSE,
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

- tlim:

  A numeric vector of length two providing lower and upper limits for
  the time axis, e.g. `c(1980, 2000)`. Use `NA` to apply no limit at
  that end. Default is `c(NA, NA)`.

- y_ticks:

  The approximate number of ticks desired on the y axis.

- ylim:

  A numeric vector of length two providing lower and upper limits for
  the y axis. Use `NA` to refer to the existing minimum or maximum. Any
  values below 1e-20 are always cut off.

- total:

  A boolean value that determines whether the total biomass from all
  species is plotted as well. Default is FALSE.

- background:

  A boolean value that determines whether background species are
  included. Ignored if the model does not contain background species.
  Default is TRUE.

- highlight:

  Name or vector of names of the species to be highlighted.

- log:

  Character string specifying which axes should use log10 scales, in the
  same form as the base
  [`plot()`](https://sizespectrum.org/mizer/reference/plot.md) argument.
  For example, `"x"`, `"y"`, `"xy"` or `""`. If supplied, this overrides
  `log_x` and `log_y`. For backward compatibility, `TRUE` and `FALSE`
  are interpreted as setting only `log_y`.

- return_data:

  A boolean value that determines whether the formatted data used for
  the plot is returned instead of the plot itself. Default is FALSE.

- log_x:

  If `TRUE`, use a log10 x-axis. Default is `FALSE`.

- log_y:

  If `TRUE`, use a log10 y-axis. Default is `TRUE`.

- use_cutoff:

  If TRUE, the `biomass_cutoff` column in the species parameters is used
  as the minimum weight for each species.

- ...:

  Arguments setting the size range over which the biomass is calculated
  (see
  [`getBiomass()`](https://sizespectrum.org/mizer/reference/getBiomass.md)):

  `min_w`

  :   Smallest weight in size range. Defaults to smallest weight in the
      model.

  `max_w`

  :   Largest weight in size range. Defaults to largest weight in the
      model.

  `min_l`

  :   Smallest length in size range. If supplied, this takes precedence
      over `min_w`.

  `max_l`

  :   Largest length in size range. If supplied, this takes precedence
      over `max_w`.

## Value

A ggplot2 object, unless `return_data = TRUE`, in which case a data
frame with the four variables 'Year', 'Biomass', 'Species', 'Legend' is
returned.

## See also

[plotting_functions](https://sizespectrum.org/mizer/reference/plotting_functions.md),
[`getBiomass()`](https://sizespectrum.org/mizer/reference/getBiomass.md)

Other plotting functions:
[`addPlot()`](https://sizespectrum.org/mizer/reference/addPlot.md),
[`animate.ArrayTimeBySpeciesBySize()`](https://sizespectrum.org/mizer/reference/animate.md),
[`plot`](https://sizespectrum.org/mizer/reference/plot.md),
[`plot2()`](https://sizespectrum.org/mizer/reference/plot2.md),
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
[`plotYieldGear()`](https://sizespectrum.org/mizer/reference/plotYieldGear.md),
[`plotting_functions`](https://sizespectrum.org/mizer/reference/plotting_functions.md)

## Examples

``` r
# \donttest{
plotBiomass(NS_sim)

plotBiomass(NS_sim, species = c("Sandeel", "Herring"), total = TRUE)

plotBiomass(NS_sim, tlim = c(1980, 1990))


# Returning the data frame
fr <- plotBiomass(NS_sim, return_data = TRUE)
str(fr)
#> 'data.frame':    528 obs. of  4 variables:
#>  $ Year   : int  1967 1968 1969 1970 1971 1972 1973 1974 1975 1976 ...
#>  $ Biomass: num  5.08e+10 5.57e+10 5.48e+10 5.32e+10 5.16e+10 ...
#>  $ Species: Factor w/ 12 levels "Sprat","Sandeel",..: 1 1 1 1 1 1 1 1 1 1 ...
#>  $ Legend : Factor w/ 12 levels "Sprat","Sandeel",..: 1 1 1 1 1 1 1 1 1 1 ...
# }
```
