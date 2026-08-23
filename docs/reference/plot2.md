# Compare two mizer arrays in a single plot

`plot2()` compares two compatible mizer array objects in a single
ggplot. Colours identify species or groups, and linetype identifies
which object the values came from.

## Usage

``` r
plot2(
  x,
  y,
  name1 = "First",
  name2 = "Second",
  species = NULL,
  log_x,
  log_y,
  log = NULL,
  ylim = c(NA, NA),
  total = FALSE,
  background = TRUE,
  highlight = NULL,
  y_ticks = 6,
  ...
)
```

## Arguments

- x:

  The first of two compatible mizer array objects to compare. Can be an
  `ArraySpeciesBySize`, `ArrayTimeBySpecies`,
  `ArrayTimeBySpeciesBySize`, `ArrayResourceBySize` or
  `ArrayTimeByResourceBySize` object.

- y:

  The second mizer array object, compatible with `x`.

- name1, name2:

  Labels for the two objects, used in the linetype legend.

- species:

  Character vector of species to include. `NULL` (default) means all
  species. A resource array holds a single spectrum, so this argument is
  not used by the resource methods, which warn if it is set.

- log_x:

  If `TRUE`, use a log10 x-axis. Default is `TRUE` for size spectra and
  `FALSE` for time series.

- log_y:

  If `TRUE`, use a log10 y-axis. Default is `FALSE` for
  `ArraySpeciesBySize` and `TRUE` for `ArrayTimeBySpecies` and for the
  resource classes.

- log:

  Character string specifying which axes should use log10 scales, in the
  same form as the base
  [`plot()`](https://sizespectrum.org/mizer/reference/plot.md) argument.
  For example, `"x"`, `"y"`, `"xy"` or `""`. If supplied, this overrides
  `log_x` and `log_y`.

- ylim:

  A numeric vector of length two providing lower and upper limits for
  the value (y) axis. Use `NA` to refer to the existing minimum or
  maximum.

- total:

  A boolean value that determines whether the total is plotted as well.
  The total is the total of everything the array holds, every species
  and every size, whatever is drawn. Default is `FALSE`. Not used by the
  resource methods, which warn if it is set.

- background:

  A boolean value that determines whether background species are
  included. Ignored if the model does not contain background species.
  Default is `TRUE`. Not used by the resource methods, which warn if it
  is set.

- highlight:

  Name or vector of names of the species to be highlighted with a
  thicker line.

- y_ticks:

  The approximate number of ticks desired on the y axis.

- ...:

  Further arguments used by only some of the methods:

  **For the `ArraySpeciesBySize`, `ArrayTimeBySpeciesBySize`,
  `ArrayResourceBySize` and `ArrayTimeByResourceBySize` methods:**

  `wlim`

  :   A numeric vector of length two providing lower and upper limits
      for the weight (x) axis. Use `NA` to refer to the existing minimum
      or maximum.

  **For the `ArraySpeciesBySize` and `ArrayTimeBySpeciesBySize`
  methods:**

  `all.sizes`

  :   If `FALSE` (default), values outside a species' size range
      (`w_min` to `w_max`) are removed.

  `llim`

  :   A numeric vector of length two providing lower and upper limits
      for the length (x) axis when `size_axis = "l"`. Use `NA` to refer
      to the existing minimum or maximum.

  `size_axis`

  :   Whether to plot size as weight (`"w"`, default) or length (`"l"`),
      using the allometric weight-length relationship of each species,
      or of the resource, see
      [`resource_params()`](https://sizespectrum.org/mizer/reference/resource_params.md).

  `per_log_size`

  :   For an array that holds a density, whether to plot it per
      logarithmic size (`TRUE`) rather than per size (`FALSE`). The
      default, `NULL`, plots the density as it stands. Unlike
      `size_axis` this needs no weight-length relationship, so it is
      available for the resource classes too. An error for an array that
      does not hold a density.

  **For `ArrayTimeBySpecies` methods:**

  `tlim`

  :   A numeric vector of length two providing lower and upper limits
      for the time axis, e.g. `c(1980, 2000)`. Use `NA` to apply no
      limit at that end. Default is `c(NA, NA)`.

  **For the `ArrayTimeBySpeciesBySize` and `ArrayTimeByResourceBySize`
  methods:**

  `time`

  :   The time to display. Default (`NULL`) is the final time step.

## Value

A ggplot2 object.

## See also

Other plotting functions:
[`addPlot()`](https://sizespectrum.org/mizer/reference/addPlot.md),
[`animate()`](https://sizespectrum.org/mizer/reference/animate.md),
[`plot`](https://sizespectrum.org/mizer/reference/plot.md),
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
[`plotYieldGear()`](https://sizespectrum.org/mizer/reference/plotYieldGear.md),
[`plotYieldVsF()`](https://sizespectrum.org/mizer/reference/plotYieldVsF.md),
[`plotting_functions`](https://sizespectrum.org/mizer/reference/plotting_functions.md)

## Examples

``` r
# \donttest{
plot2(getEncounter(NS_params), getEncounter(NS_params))

plot2(getResourceMort(NS_params), getResourceMort(NS_params))

# }
```
