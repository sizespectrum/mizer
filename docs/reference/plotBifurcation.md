# Draw a bifurcation diagram over fishing effort

**\[experimental\]**

Sweeps the fishing effort over a range of values and, for each value,
runs the full dynamics to their attractor and records the long-term
range of a summary quantity (biomass by default). The result is a
bifurcation diagram with fishing effort on the x-axis.

For each effort value the attractor is found in two stages. First
[`projectToSteady()`](https://sizespectrum.org/mizer/reference/projectToSteady.md)
runs the dynamics until they settle, stopping early once it detects a
stable steady state or a limit cycle and reporting which via its
`"convergence"` attribute. Then the settled state is projected forward
once more with
[`project()`](https://sizespectrum.org/mizer/reference/project.md) over
a short sampling window — one full period for a limit cycle, or a few
years otherwise — and the minimum and maximum of the chosen quantity
over that window are taken as the attractor envelope.

For a stable steady state the minimum and maximum coincide and the
species is drawn as a single line; once the dynamics settle onto a limit
cycle the two separate and the species is drawn as a shaded band between
them. The onset of the band therefore marks the Hopf bifurcation (see
[`vignette("dynamic_stability")`](https://sizespectrum.org/mizer/articles/dynamic_stability.md)).

By default the sweep uses *continuation*: the projection at each effort
value starts from the attractor reached at the previous value rather
than from the original state. This shortens the transient and keeps the
sweep on a single branch of the attractor. Because it follows one
branch, the diagram can look different for an increasing versus a
decreasing effort sequence if the model has coexisting attractors
(hysteresis); pass a decreasing `effort` vector to trace the other
direction.

## Usage

``` r
plotBifurcation(
  params,
  effort = seq(0, 2, length.out = 21),
  species = NULL,
  value = c("biomass", "yield", "ssb"),
  t_max = 100,
  t_sample = NULL,
  t_sample_default = 10,
  t_save = 0.25,
  tol = 0.01,
  amplitude_tol = 0.01,
  extinction_threshold = 1e-06,
  continuation = TRUE,
  return_data = FALSE,
  progress_bar = TRUE,
  ytrans = "log10"
)
```

## Arguments

- params:

  A
  [MizerParams](https://sizespectrum.org/mizer/reference/MizerParams-class.md)
  object.

- effort:

  A numeric vector of fishing effort values for the x-axis. The same
  effort is applied to every gear (as in `project(effort = value)`).
  Default `seq(0, 2, length.out = 21)`.

- species:

  The species to include. By default all target species. A vector of
  species names or indices, as for other mizer plotting functions.

- value:

  The quantity for the y-axis, one of `"biomass"` (default), `"yield"`
  or `"ssb"`, computed with
  [`getBiomass()`](https://sizespectrum.org/mizer/reference/getBiomass.md),
  [`getYield()`](https://sizespectrum.org/mizer/reference/getYield.md)
  or [`getSSB()`](https://sizespectrum.org/mizer/reference/getSSB.md)
  respectively.

- t_max:

  The maximum number of years to run the settling stage
  ([`projectToSteady()`](https://sizespectrum.org/mizer/reference/projectToSteady.md))
  at each effort value.

- t_sample:

  The length in years of the window over which the settled attractor is
  sampled to measure the envelope. If `NULL` (default) it is chosen
  automatically: one full period for a detected limit cycle, or
  `t_sample_default` years otherwise.

- t_sample_default:

  The sampling window used when no limit-cycle period is available (a
  stable or non-converged run). Default `10`.

- t_save:

  The interval at which the sampling window is saved, controlling how
  finely the cycle envelope is resolved. Default `0.25`.

- tol:

  Convergence tolerance for the settling stage, passed to
  [`projectToSteady()`](https://sizespectrum.org/mizer/reference/projectToSteady.md).
  Tighter (smaller) values settle the stable branch more fully, giving
  cleaner single lines below the bifurcation at the cost of longer runs.
  Default `0.01`.

- amplitude_tol:

  The minimum relative biomass amplitude for a run to be classified as a
  limit cycle, passed to
  [`projectToSteady()`](https://sizespectrum.org/mizer/reference/projectToSteady.md).
  Default `0.01`.

- extinction_threshold:

  The relative reproduction collapse below which a species is treated as
  extinct, passed to
  [`projectToSteady()`](https://sizespectrum.org/mizer/reference/projectToSteady.md).
  Default `1e-6`.

- continuation:

  If `TRUE` (default) each settling run warm-starts from the attractor
  of the previous effort value.

- return_data:

  If `TRUE` the data frame underlying the plot is returned instead of
  the plot. Default `FALSE`.

- progress_bar:

  If `TRUE` (default) a text progress bar is shown while the effort
  values are swept.

- ytrans:

  Transformation for the y-axis, `"log10"` (default) or `"identity"`.

## Value

A ggplot2 object, or a data frame with columns `Effort`, `Species`,
`ymin`, `ymax` and `type` (the attractor type reported by
[`projectToSteady()`](https://sizespectrum.org/mizer/reference/projectToSteady.md))
if `return_data = TRUE`.

## See also

[`getStability()`](https://sizespectrum.org/mizer/reference/getStability.md),
[`projectToSteady()`](https://sizespectrum.org/mizer/reference/projectToSteady.md),
[`plotBiomass()`](https://sizespectrum.org/mizer/reference/plotBiomass.md)

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

## Examples

``` r
# \donttest{
plotBifurcation(NS_params, effort = seq(0, 2, length.out = 11))
#>   |                                                                              |                                                                      |   0%  |                                                                              |======                                                                |   9%  |                                                                              |=============                                                         |  18%  |                                                                              |===================                                                   |  27%  |                                                                              |=========================                                             |  36%  |                                                                              |================================                                      |  45%  |                                                                              |======================================                                |  55%  |                                                                              |=============================================                         |  64%  |                                                                              |===================================================                   |  73%  |                                                                              |=========================================================             |  82%  |                                                                              |================================================================      |  91%  |                                                                              |======================================================================| 100%

# }
```
