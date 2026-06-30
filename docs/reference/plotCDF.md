# Plot cumulative abundance or biomass distributions

`plotCDF()` plots the cumulative distribution over body size from small
to large sizes. It uses the same spectra data preparation as
[`plotSpectra()`](https://sizespectrum.org/mizer/reference/plotSpectra.md).
The density is first multiplied by `w^power`, then integrated over size.
With `normalise = TRUE`, each curve is divided by its final value so
that it ends at 1.

## Usage

``` r
plotCDF(
  object,
  species = NULL,
  wlim = c(NA, NA),
  llim = c(NA, NA),
  ylim = c(NA, NA),
  power = 1,
  biomass = TRUE,
  total = FALSE,
  resource = FALSE,
  background = TRUE,
  highlight = NULL,
  normalise = TRUE,
  log_x = TRUE,
  log_y = FALSE,
  log = NULL,
  size_axis = c("w", "l"),
  return_data = FALSE,
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
  is FALSE.

- background:

  A boolean value that determines whether background species are
  included. Ignored if the model does not contain background species.
  Default is TRUE.

- highlight:

  Name or vector of names of the species to be highlighted by being
  plotted with thicker lines.

- normalise:

  If `TRUE` (default), plot the cumulative proportion. If `FALSE`, plot
  the cumulative abundance, biomass, or other unnormalised integral.

- log_x:

  If `TRUE` (default), use a log10 x-axis.

- log_y:

  If `TRUE`, use a log10 y-axis. Default is `FALSE`.

- log:

  Character string specifying which axes should use a log10 scale, in
  the same form as the base
  [`plot()`](https://sizespectrum.org/mizer/reference/plot.md) argument.
  If supplied, this overrides `log_x` and `log_y`.

- size_axis:

  Whether to plot size as weight (`"w"`, default) or length (`"l"`),
  using the allometric weight-length relationship.

- return_data:

  A boolean value that determines whether the formatted data used for
  the plot is returned instead of the plot itself. Default is FALSE.

- ...:

  Further arguments used by only some of the methods:

  **For `MizerSim` methods:**

  - `time_range`: The time range (either a vector of values, a vector of
    min and max time, or a single value) to average the abundances over.
    Default is the final time step.

  - `geometric_mean`: **\[experimental\]** If `TRUE` then the average of
    the abundances over the time range is a geometric mean instead of
    the default arithmetic mean.

## Value

A ggplot2 object, unless `return_data = TRUE`, in which case a data
frame with the four variables 'w' (or 'l' if `size_axis = "l"`),
'value', 'Species', 'Legend' is returned. `plotlyCDF()` returns a plotly
object.

## Details

`plotlyCDF()` is the interactive plotly version. To compare cumulative
distributions from two objects, use
[`plotCDF2()`](https://sizespectrum.org/mizer/reference/plotCDF2.md).

## See also

[`plotSpectra()`](https://sizespectrum.org/mizer/reference/plotSpectra.md),
[`plotCDF2()`](https://sizespectrum.org/mizer/reference/plotCDF2.md)

Other plotting functions:
[`addPlot()`](https://sizespectrum.org/mizer/reference/addPlot.md),
[`animate.ArrayTimeBySpeciesBySize()`](https://sizespectrum.org/mizer/reference/animate.md),
[`plot`](https://sizespectrum.org/mizer/reference/plot.md),
[`plot2()`](https://sizespectrum.org/mizer/reference/plot2.md),
[`plotBiomass()`](https://sizespectrum.org/mizer/reference/plotBiomass.md),
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
plotCDF(NS_params, species = c("Cod", "Herring"))

plotCDF(NS_sim, power = 0, normalise = FALSE)

# }
```
