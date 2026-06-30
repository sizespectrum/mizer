# Plot diet, resolved by prey species, as function of predator at size.

**\[experimental\]** Plots the proportions with which each prey species
contributes to the total biomass consumed by the specified predator
species, as a function of the predator's size. These proportions are
obtained with
[`getDiet()`](https://sizespectrum.org/mizer/reference/getDiet.md).

## Usage

``` r
plotDiet(
  object,
  species = NULL,
  wlim = c(NA, NA),
  llim = c(NA, NA),
  size_axis = c("w", "l"),
  return_data = FALSE,
  log_x = TRUE,
  log_y = FALSE,
  log = NULL,
  ...
)
```

## Arguments

- object:

  An object of class
  [MizerSim](https://sizespectrum.org/mizer/reference/MizerSim-class.md)
  or
  [MizerParams](https://sizespectrum.org/mizer/reference/MizerParams-class.md).

- species:

  The species to be selected. Optional. By default all target species
  are selected. A vector of species names, or a numeric vector with the
  species indices, or a logical vector indicating for each species
  whether it is to be selected (TRUE) or not.

- wlim:

  A numeric vector of length two providing lower and upper limits for
  the weight (x) axis. Use `NA` to auto-scale to the data range.

- llim:

  A numeric vector of length two providing lower and upper limits for
  the length (x) axis when `size_axis = "l"`. Use `NA` to auto-scale to
  the data range.

- size_axis:

  Whether to plot size as weight (`"w"`, default) or length (`"l"`),
  using the allometric weight-length relationship.

- return_data:

  A boolean value that determines whether the formatted data used for
  the plot is returned instead of the plot itself. Default is FALSE.

- log_x:

  If `TRUE` (default), use a log10 x-axis.

- log_y:

  If `TRUE`, use a log10 y-axis. Default is `FALSE`.

- log:

  Character string specifying which axes should use log10 scales, in the
  same form as the base
  [`plot()`](https://sizespectrum.org/mizer/reference/plot.md) argument.
  For example, `"x"`, `"y"`, `"xy"` or `""`. If supplied, this overrides
  `log_x` and `log_y`.

- ...:

  Further arguments used by only some of the methods:

  **For `MizerSim` methods:**

  - `time_range`: The time range (either a vector of values, a vector of
    min and max time, or a single value) over which to average the diet.
    The consumption rates are averaged over this range and then
    normalised to proportions. Default is the final time step.

## Value

A ggplot2 object, unless `return_data = TRUE`, in which case a data
frame with the four variables 'Predator', 'w' (or 'l' if
`size_axis = "l"`), 'Proportion', 'Prey' is returned.

`plotlyDiet()` returns a plotly object.

## Details

Prey species that contribute less than 1 permille to the diet are
suppressed in the plot. The plot only extends to predator sizes where
the predator has a meaningful abundance (defined as having a biomass
density greater than 0.1% of its maximum biomass density).

If more than one predator species is selected, then the plot contains
one facet for each species.

## See also

[`getDiet()`](https://sizespectrum.org/mizer/reference/getDiet.md)

Other plotting functions:
[`addPlot()`](https://sizespectrum.org/mizer/reference/addPlot.md),
[`animate.ArrayTimeBySpeciesBySize()`](https://sizespectrum.org/mizer/reference/animate.md),
[`plot`](https://sizespectrum.org/mizer/reference/plot.md),
[`plot2()`](https://sizespectrum.org/mizer/reference/plot2.md),
[`plotBiomass()`](https://sizespectrum.org/mizer/reference/plotBiomass.md),
[`plotCDF()`](https://sizespectrum.org/mizer/reference/plotCDF.md),
[`plotCDF2()`](https://sizespectrum.org/mizer/reference/plotCDF2.md),
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
plotDiet(NS_params, species = "Cod")

plotDiet(NS_params, species = 5:9)


# Returning the data frame
fr <- plotDiet(NS_params, species = "Cod", return_data = TRUE)
str(fr)
#> 'data.frame':    800 obs. of  4 variables:
#>  $ Predator  : Factor w/ 1 level "Cod": 1 1 1 1 1 1 1 1 1 1 ...
#>  $ w         : num  0.486 0.58 0.693 0.827 0.987 1.18 1.4 1.68 2 2.39 ...
#>  $ Prey      : Factor w/ 14 levels "External","Resource",..: 14 14 14 14 14 14 14 14 14 14 ...
#>  $ Proportion: num  0.00113 0.00139 0.0017 0.00208 0.00254 ...
# }
```
