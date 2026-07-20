# Plot mizer arrays

Many mizer functions return values that depend on species and either
size or time. `plot()` creates a ggplot2 figure with one line for each
species showing the values against size or against time (depending on
the type of output).
[`plotHover()`](https://sizespectrum.org/mizer/reference/plotHover.md)
creates an interactive version of the same figure.

## Details

This works because the mizer functions that give values that depend on
species and size return an `ArraySpeciesBySize` object and those that
give values that depend on species and time return an
`ArrayTimeBySpecies` object. These objects have attributes that store
the name of the value, its units, and a reference to the `MizerParams`
object that the value was computed from. This allows the plots to be
automatically labelled and coloured appropriately.

To compare two mizer arrays in a single plot, use
[`plot2()`](https://sizespectrum.org/mizer/reference/plot2.md). To show
the relative difference between two arrays, use
[`plotRelative()`](https://sizespectrum.org/mizer/reference/plotRelative.md).

All methods return a ggplot2 object, unless `return_data = TRUE`, in
which case they return the underlying data frame instead.
[`plotHover()`](https://sizespectrum.org/mizer/reference/plotHover.md)
returns a plotly object.

Arguments used by all methods:

- `species`:

  Character vector of species to include. `NULL` (default) means all
  species.

- `highlight`:

  Name or vector of names of the species to be highlighted.

- `total`:

  A boolean value that determines whether the total over all selected
  species is plotted as well. Default is `FALSE`.

- `background`:

  A boolean value that determines whether background species are
  included. Ignored if the model does not contain background species.
  Default is `TRUE`.

- `return_data`:

  If `TRUE`, return the data frame instead of the plot.

- `log_x`:

  If `TRUE`, use a log10 x-axis. The default depends on the method; see
  its own help page.

- `log_y`:

  If `TRUE`, use a log10 y-axis. The default depends on the method; see
  its own help page.

- `log`:

  Character string specifying which axes should use log10 scales, in the
  same form as the base `plot()` argument. For example, `"x"`, `"y"`,
  `"xy"` or `""`. If supplied, this overrides `log_x` and `log_y`.

- `ylim`:

  A numeric vector of length two providing lower and upper limits for
  the value (y) axis. Use `NA` to refer to the existing minimum or
  maximum.

- `y_ticks`:

  The approximate number of ticks desired on the y axis.

Additional arguments for
[`plot.ArraySpeciesBySize()`](https://sizespectrum.org/mizer/reference/plot.ArraySpeciesBySize.md)
and
[`plot.ArrayTimeBySpeciesBySize()`](https://sizespectrum.org/mizer/reference/plot.ArrayTimeBySpeciesBySize.md):

- `all.sizes`:

  If `FALSE` (default), values outside a species' size range (`w_min` to
  `w_max`) are removed.

- `wlim`:

  A numeric vector of length two providing lower and upper limits for
  the weight (x) axis. Use `NA` to refer to the existing minimum or
  maximum.

- `llim`:

  A numeric vector of length two providing lower and upper limits for
  the length (x) axis when `size_axis = "l"`. Use `NA` to refer to the
  existing minimum or maximum.

- `size_axis`:

  Whether to plot size as weight (`"w"`, default) or length (`"l"`),
  using the allometric weight-length relationship.

Additional argument for
[`plot.ArrayTimeBySpecies()`](https://sizespectrum.org/mizer/reference/plot.ArrayTimeBySpecies.md):

- `tlim`:

  A numeric vector of length two providing lower and upper limits for
  the time axis, e.g. `c(1980, 2000)`. Use `NA` to apply no limit at
  that end. Default is `c(NA, NA)`.

Additional argument for
[`plot.ArrayTimeBySpeciesBySize()`](https://sizespectrum.org/mizer/reference/plot.ArrayTimeBySpeciesBySize.md)
and
[`plot.ArrayTimeByResourceBySize()`](https://sizespectrum.org/mizer/reference/plot.ArrayTimeByResourceBySize.md):

- `time`:

  The time to display. Default (`NULL`) is the final time step.

See the individual method help pages for each method's exact arguments
and defaults:
[`plot.ArraySpeciesBySize()`](https://sizespectrum.org/mizer/reference/plot.ArraySpeciesBySize.md),
[`plot.ArrayTimeBySpecies()`](https://sizespectrum.org/mizer/reference/plot.ArrayTimeBySpecies.md),
[`plot.ArrayTimeBySpeciesBySize()`](https://sizespectrum.org/mizer/reference/plot.ArrayTimeBySpeciesBySize.md),
[`plot.ArrayResourceBySize()`](https://sizespectrum.org/mizer/reference/plot.ArrayResourceBySize.md),
[`plot.ArrayTimeByResourceBySize()`](https://sizespectrum.org/mizer/reference/plot.ArrayTimeByResourceBySize.md).

## See also

Other plotting functions:
[`addPlot()`](https://sizespectrum.org/mizer/reference/addPlot.md),
[`animate.ArrayTimeBySpeciesBySize()`](https://sizespectrum.org/mizer/reference/animate.md),
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
plot(getEncounter(NS_params))

plot(getFeedingLevel(NS_params), species = c("Cod", "Herring"))

plot(getPredMort(NS_params), species = c("Cod", "Herring"),
     size_axis = "l")

plot(getBiomass(NS_sim))

plot(getBiomass(NS_sim), species = c("Cod", "Herring"), total = TRUE)

plot(getYield(NS_sim), species = c("Cod", "Herring"))

plot(getFMort(NS_sim), time = 2010)

plot(getResourceMort(NS_params))

plot(initialNResource(NS_params))
#> Warning: log-10 transformation introduced infinite values.

plot(NResource(NS_sim))
#> Warning: log-10 transformation introduced infinite values.

# }
```
