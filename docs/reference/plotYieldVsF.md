# Plot the yield of a species against the fishing mortality on it

**\[experimental\]**

Varies the fishing mortality on one species over a range of values,
leaving the fishing on every other species unchanged, and plots the
long-term yield of that species against it. The fishing mortality at
which the yield is largest is \\F\_{MSY}\\, and is marked on the plot by
default.

This is
[`scanModel()`](https://sizespectrum.org/mizer/reference/scanModel.md)
with
[`scanFishingMortality()`](https://sizespectrum.org/mizer/reference/scanEffort.md)
as its setter and
[`getYield()`](https://sizespectrum.org/mizer/reference/getYield.md) as
the quantity it measures. Use
[`scanModel()`](https://sizespectrum.org/mizer/reference/scanModel.md)
directly to vary something other than the fishing mortality on a single
species, to measure something other than the yield, or to follow more
than one species at once.

At each fishing mortality the model is projected until it settles, and
what is plotted depends on what it settled on. At a fixed point the
yield is read straight off the settled state. On a limit cycle it is
averaged over exactly one period, and the band around the line shows the
range the yield covers over that cycle, so an oscillation is displayed
rather than silently averaged away. Fishing mortalities at which the
model settled on neither are marked with a cross and should not be
relied on; raise `t_max` for those.

The scan starts from the fishing mortality the model currently sits at
and works outwards in both directions, each arm warm-starting from the
attractor reached at the previous value.

## Usage

``` r
plotYieldVsF(
  params,
  species,
  F_range,
  F_min = 0,
  F_max = 1.5,
  no_steps = 16,
  gear = NULL,
  style = "ribbon",
  mark_max = TRUE,
  reference_lines = TRUE,
  log_y = FALSE,
  log = NULL,
  return_data = FALSE,
  progress_bar = interactive(),
  ...
)
```

## Arguments

- params:

  A
  [MizerParams](https://sizespectrum.org/mizer/reference/MizerParams-class.md)
  object.

- species:

  The name of the species whose fishing mortality is varied. Only one
  species at a time.

- F_range:

  A numeric vector of fishing mortalities for the x-axis. If missing it
  is built as `seq(F_min, F_max, length.out = no_steps)`.

- F_min, F_max, no_steps:

  Used to build `F_range` when that is missing.

- gear:

  The name of the gear whose fishing mortality on the species is varied.
  Only needed when several gears catch the species; if NULL (default),
  the fishing mortality from all of them is replaced. See
  [`scanFishingMortality()`](https://sizespectrum.org/mizer/reference/scanEffort.md).

- style:

  How the range covered on a limit cycle is drawn, see
  [`plot.MizerScan()`](https://sizespectrum.org/mizer/reference/plot.MizerScan.md).
  The default `"ribbon"` draws the average as a line inside the band.

- mark_max:

  Whether to mark the fishing mortality at which the yield is largest,
  which is \\F\_{MSY}\\. Default TRUE.

- reference_lines:

  Whether to draw the `F_MSY` species parameter, if the species has one,
  as a vertical line. See
  [`plot.MizerScan()`](https://sizespectrum.org/mizer/reference/plot.MizerScan.md).

- log_y, log:

  Whether to use a logarithmic y-axis, see
  [`parsePlotLog()`](https://sizespectrum.org/mizer/reference/parsePlotLog.md).
  Unlike most mizer plots this defaults to FALSE, because the yield is
  exactly zero at zero fishing mortality and a logarithmic axis would
  drop the point that anchors the curve.

- return_data:

  If TRUE the
  [MizerScan](https://sizespectrum.org/mizer/reference/MizerScan.md)
  object underlying the plot is returned instead of the plot. Default
  FALSE.

- progress_bar:

  If TRUE a text progress bar is shown while the fishing mortalities are
  swept. Defaults to
  [`interactive()`](https://rdrr.io/r/base/interactive.html).

- ...:

  Further arguments are passed on to
  [`scanModel()`](https://sizespectrum.org/mizer/reference/scanModel.md).

## Value

A ggplot2 object, or, if `return_data = TRUE`, the
[MizerScan](https://sizespectrum.org/mizer/reference/MizerScan.md)
object holding the data. The fishing mortality giving the largest yield
is available from that object as `attr(scan, "at_max")`.

## See also

[`scanModel()`](https://sizespectrum.org/mizer/reference/scanModel.md),
[`scanFishingMortality()`](https://sizespectrum.org/mizer/reference/scanEffort.md),
[`plot.MizerScan()`](https://sizespectrum.org/mizer/reference/plot.MizerScan.md),
[`getYield()`](https://sizespectrum.org/mizer/reference/getYield.md)

Other plotting functions:
[`addPlot()`](https://sizespectrum.org/mizer/reference/addPlot.md),
[`animate()`](https://sizespectrum.org/mizer/reference/animate.md),
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

Other scan functions:
[`MizerScan()`](https://sizespectrum.org/mizer/reference/MizerScan.md),
[`plot.MizerScan()`](https://sizespectrum.org/mizer/reference/plot.MizerScan.md),
[`scanEffort()`](https://sizespectrum.org/mizer/reference/scanEffort.md),
[`scanModel()`](https://sizespectrum.org/mizer/reference/scanModel.md)

## Examples

``` r
# \donttest{
plotYieldVsF(NS_params, "Cod", F_max = 1.5, no_steps = 8)


# The fishing mortality that maximises the yield
scan <- plotYieldVsF(NS_params, "Cod", F_max = 1.5, no_steps = 8,
                     return_data = TRUE)
attr(scan, "at_max")
#>       Cod 
#> 0.8571429 
# }
```
