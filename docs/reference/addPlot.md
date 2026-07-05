# Add lines to an existing plot

**\[experimental\]** `addPlot()` adds another set of values to an
existing ggplot. The first method supports adding an
`ArraySpeciesBySize` object to a compatible plot, for example to compare
the same rate before and after a model change. The method checks whether
the existing plot uses a compatible x variable, and warns if the y
variable or y-axis units appear to differ.

## Usage

``` r
addPlot(
  plot,
  x,
  species = NULL,
  total = FALSE,
  background = TRUE,
  colour = NULL,
  linetype = "dashed",
  linewidth = 0.8,
  alpha = 1,
  ...
)
```

## Arguments

- plot:

  A ggplot2 object to which the new values should be added.

- x:

  An object containing the values to add.

- species:

  Character vector of species to include. `NULL` (default) means all
  species.

- total:

  A boolean value that determines whether the total over all selected
  species is plotted as well. Default is `FALSE`.

- background:

  A boolean value that determines whether background species are
  included. Ignored if the model does not contain background species.
  Default is `TRUE`.

- colour:

  Optional fixed colour for the added lines. If `NULL`, the species
  colours from the existing plot are used.

- linetype:

  Optional fixed line type for the added lines. If `NULL`, the species
  line types from the existing plot are used.

- linewidth:

  Width of the added lines.

- alpha:

  Transparency of the added lines.

- ...:

  Further arguments used by only some of the methods:

  **For `ArraySpeciesBySize` methods:**

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

  `ylim`

  :   A numeric vector of length two providing lower and upper limits
      for the value (y) axis.

## Value

A ggplot2 object.

## See also

Other plotting functions:
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
[`plotYieldGear()`](https://sizespectrum.org/mizer/reference/plotYieldGear.md),
[`plotting_functions`](https://sizespectrum.org/mizer/reference/plotting_functions.md)

## Examples

``` r
# \donttest{
p <- plot(getEncounter(NS_params), species = "Cod")
addPlot(p, getEncounter(NS_params), species = "Cod")

# }
```
