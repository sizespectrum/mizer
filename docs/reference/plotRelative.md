# Plot relative difference between two mizer arrays

`plotRelative()` plots the difference between two compatible mizer array
objects relative to their average. If the values in the first object are
\\N_1\\ and the values in the second are \\N_2\\, it plots \$\$2 (N_2 -
N_1) / (N_1 + N_2).\$\$

## Usage

``` r
plotRelative(
  x,
  y,
  species = NULL,
  log_x,
  ylim = c(NA, NA),
  total = FALSE,
  background = TRUE,
  ...
)
```

## Arguments

- x:

  The first of two compatible mizer array objects to compare. Can be an
  `ArraySpeciesBySize`, `ArrayTimeBySpecies`, or
  `ArrayTimeBySpeciesBySize` object.

- y:

  The second mizer array object, compatible with `x`.

- species:

  Character vector of species to include. `NULL` (default) means all
  species.

- log_x:

  If `TRUE`, use a log10 x-axis. Default is `TRUE` for size spectra and
  `FALSE` for time series.

- ylim:

  A numeric vector of length two providing lower and upper limits for
  the value (y) axis.

- total:

  A boolean value that determines whether the total over all selected
  species is plotted as well. Default is `FALSE`.

- background:

  A boolean value that determines whether background species are
  included. Ignored if the model does not contain background species.
  Default is `TRUE`.

- ...:

  Further arguments used by only some of the methods:

  **For `ArraySpeciesBySize` and `ArrayTimeBySpeciesBySize` methods:**

  `all.sizes`

  :   If `FALSE` (default), values outside a species' size range
      (`w_min` to `w_max`) are removed.

  `wlim`

  :   A numeric vector of length two providing lower and upper limits
      for the weight (x) axis. Use `NA` to refer to the existing minimum
      or maximum.

  `llim`

  :   A numeric vector of length two providing lower and upper limits
      for the length (x) axis when `size_axis = "l"`. Use `NA` to refer
      to the existing minimum or maximum.

  `size_axis`

  :   Whether to plot size as weight (`"w"`, default) or length (`"l"`),
      using the allometric weight-length relationship.

  **For `ArrayTimeBySpecies` methods:**

  `tlim`

  :   A numeric vector of length two providing lower and upper limits
      for the time axis, e.g. `c(1980, 2000)`. Use `NA` to apply no
      limit at that end. Default is `c(NA, NA)`.

  **For `ArrayTimeBySpeciesBySize` methods:**

  `time`

  :   The time to display. Default (`NULL`) is the final time step.

## Value

A ggplot2 object.

## See also

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
[`plotSpectra()`](https://sizespectrum.org/mizer/reference/plotSpectra.md),
[`plotSpectra2()`](https://sizespectrum.org/mizer/reference/plotSpectra2.md),
[`plotSpectraRelative()`](https://sizespectrum.org/mizer/reference/plotSpectraRelative.md),
[`plotYield()`](https://sizespectrum.org/mizer/reference/plotYield.md),
[`plotYieldGear()`](https://sizespectrum.org/mizer/reference/plotYieldGear.md),
[`plotting_functions`](https://sizespectrum.org/mizer/reference/plotting_functions.md)

## Examples

``` r
# \donttest{
params <- NS_params
given_species_params(params)["Cod", "w_mat"] <- 1200
plotRelative(getEGrowth(NS_params), getEGrowth(params),
             wlim = c(500, 2000), log_x = FALSE, species = "Cod")

# }
```
