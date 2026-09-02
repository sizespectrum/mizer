# Guide: Analysing and plotting mizer results

This guide gives an overview of the functions available in mizer for
analysing the results of simulations and creating plots. For full
documentation of each function, follow the links.

Mizer ships a large family of extraction, summary, and plotting
functions. **Always prefer these over hand-written array wrangling or
custom ggplot code** — they handle size-range integration, species
colours/linetypes, and units for you.

Most functions accept **either** a
[`MizerSim`](https://sizespectrum.org/mizer/reference/MizerSim.md)
object (returning a **time series**) **or** a
[`MizerParams`](https://sizespectrum.org/mizer/reference/MizerParams.md)
object (returning a **single value** from the initial state). So
[`getBiomass(sim)`](https://sizespectrum.org/mizer/reference/getBiomass.md)
gives biomass over time, `getBiomass(params)` gives biomass now.

To get the single value **at one time step** of a simulation, extract a
`MizerParams` snapshot with
[`finalParams(sim)`](https://sizespectrum.org/mizer/reference/getParams.md)
(last step),
[`initialParams(sim)`](https://sizespectrum.org/mizer/reference/getParams.md)
(first step), or
[`getParams(sim, time_range = ...)`](https://sizespectrum.org/mizer/reference/getParams.md)
(averaged over a range) and pass that in:

``` r

getMeanMaxWeight(finalParams(sim))             # value at the last time step
getSSB(getParams(sim, time_range = 1990:2000)) # averaged over a period
```

------------------------------------------------------------------------

## Accessing simulation arrays

These extract raw arrays from a `MizerSim` object.

| Function | Returns | Dimensions |
|----|----|----|
| [`N(sim)`](https://sizespectrum.org/mizer/reference/N.md) | species abundance density | time × species × size |
| [`NResource(sim)`](https://sizespectrum.org/mizer/reference/N.md) | resource abundance density | time × size |
| [`finalN(sim)`](https://sizespectrum.org/mizer/reference/finalN.md) | species abundance at last time | species × size |
| [`finalNResource(sim)`](https://sizespectrum.org/mizer/reference/finalN.md) | resource abundance at last time | size |
| [`getEffort(sim)`](https://sizespectrum.org/mizer/reference/getEffort.md) | fishing effort | time × gear |
| [`getTimes(sim)`](https://sizespectrum.org/mizer/reference/getTimes.md) | saved time steps | time |

``` r

N(sim)[, , 1]              # time × species in smallest size class
N(sim)["2010", "Cod", ]    # size vector for Cod in year 2010
finalN(sim)["Cod", ]       # size vector for Cod at the final time step
```

------------------------------------------------------------------------

## Summary functions

These functions compute derived quantities from abundances. All accept
`MizerSim` or `MizerParams`. The result is a classed array that can be
plotted directly with
[`plot()`](https://sizespectrum.org/mizer/reference/plot.md) — see
below.

| Function | Returns | Dimensions |
|----|----|----|
| [`getBiomass(sim, min_w, max_w)`](https://sizespectrum.org/mizer/reference/getBiomass.md) | total biomass | time × species |
| [`getSSB(sim)`](https://sizespectrum.org/mizer/reference/getSSB.md) | spawning stock biomass | time × species |
| [`getN(sim, min_w, max_w)`](https://sizespectrum.org/mizer/reference/getN.md) | total abundance | time × species |
| [`getYield(sim)`](https://sizespectrum.org/mizer/reference/getYield.md) | total yield across gears | time × species |
| [`getYieldGear(sim)`](https://sizespectrum.org/mizer/reference/getYieldGear.md) | yield by gear | time × gear × species |
| [`getFeedingLevel(sim)`](https://sizespectrum.org/mizer/reference/getFeedingLevel.md) | feeding level at size | time × species × size |
| [`getPredMort(sim)`](https://sizespectrum.org/mizer/reference/getPredMort.md) | predation mortality at size | time × species × size |
| [`getFMort(sim)`](https://sizespectrum.org/mizer/reference/getFMort.md) | fishing mortality at size | time × species × size |
| [`getFMortGear(sim)`](https://sizespectrum.org/mizer/reference/getFMortGear.md) | fishing mortality by gear | time × gear × species × size |
| [`getDiet(params)`](https://sizespectrum.org/mizer/reference/getDiet.md) | diet resolved by prey at size | predator × size × prey |
| [`getTrophicLevel(params)`](https://sizespectrum.org/mizer/reference/getTrophicLevel.md) | trophic level at size | species × size |
| [`getTrophicLevelBySpecies(params)`](https://sizespectrum.org/mizer/reference/getTrophicLevelBySpecies.md) | mean trophic level per species | species |

**Size range:**
[`getBiomass()`](https://sizespectrum.org/mizer/reference/getBiomass.md)
and [`getN()`](https://sizespectrum.org/mizer/reference/getN.md) accept
`min_w`, `max_w`, `min_l`, `max_l` to restrict the calculation to a size
range.

``` r

getSSB(sim)                              # SSB of all species over time
getBiomass(sim, min_w = 10, max_w = 1e4) # biomass of 10g–10kg fish
getYield(sim)["2010", ]                  # yield in year 2010
```

------------------------------------------------------------------------

## Indicator functions

These compute community-level indicators. All accept `MizerSim` (time
series) or `MizerParams` (single value from the initial state). See
[`?indicator_functions`](https://sizespectrum.org/mizer/reference/indicator_functions.md).

| Function | Key arguments | Returns |
|----|----|----|
| [`getProportionOfLargeFish(sim)`](https://sizespectrum.org/mizer/reference/getProportionOfLargeFish.md) | `threshold_w = 100`, `biomass_proportion` | proportion of large fish through time |
| [`getMeanWeight(sim)`](https://sizespectrum.org/mizer/reference/getMeanWeight.md) | `min_w`, `max_w`, `species` | mean community weight through time |
| [`getMeanLength(sim)`](https://sizespectrum.org/mizer/reference/getMeanWeight.md) | `min_w`, `max_w`, `species` | mean community length through time |
| [`getMeanMaxWeight(sim)`](https://sizespectrum.org/mizer/reference/getMeanMaxWeight.md) | `measure = "both"/"numbers"/"biomass"` | mean asymptotic weight through time |
| [`getCommunitySlope(sim)`](https://sizespectrum.org/mizer/reference/getCommunitySlope.md) | `min_w`, `max_w`, `species` | slope, intercept, R² through time |

``` r

lfi <- getProportionOfLargeFish(sim, min_w = 10, max_w = 5000, threshold_w = 500)
slope <- getCommunitySlope(sim, min_w = 10, max_w = 5000)
```

------------------------------------------------------------------------

## Writing your own indicator

First check that a built-in does not already cover it: most custom
indicators turn out to be
[`getBiomass()`](https://sizespectrum.org/mizer/reference/getBiomass.md)/[`getN()`](https://sizespectrum.org/mizer/reference/getN.md)
over a size range, or one of the four above with different arguments. If
none fits, an indicator is an integral over the size spectrum, \\\int
N_i(w)\\ K_i(w)\\ dw\\, where \\K_i(w)\\ is a **weighting factor**
(supplied to the `weighting` argument of
[`sizeIntegral()`](https://sizespectrum.org/mizer/reference/sizeIntegral.md)).
**[`sizeIntegral()`](https://sizespectrum.org/mizer/reference/sizeIntegral.md)
does that integral for you**:

``` r

# Abundance between 10g and 5kg (default weighting factor weighting = 1)
sizeIntegral(params, min_w = 10, max_w = 5000)

# Biomass between 10g and 5kg (weighting factor is body weight)
sizeIntegral(params, weighting = w(params), min_w = 10, max_w = 5000)

# Biomass through time, wrapped ready to plot
sizeIntegral(sim, weighting = w(params), value_name = "Biomass", units = "g")
```

Give it the object and the weighting factor \\K\\ (e.g. body weight
[`w(params)`](https://sizespectrum.org/mizer/reference/w.md) for
biomass, or 1 for numbers); it selects the size range, uses the
quadrature scheme the model is actually on and wraps the result in the
appropriate mizer array class. Doing the sum by hand instead means
getting all of that right yourself, silently and only for some users.
Three things are worth knowing:

- **The size range is an argument, not a subsetting step.**
  `min_w`/`max_w` or `min_l`/`max_l` are passed through to
  [`get_size_range_array()`](https://sizespectrum.org/mizer/reference/get_size_range_array.md),
  which does the length-weight conversion per species and accepts either
  a single number or one value per species. Never subset the size grid
  by hand.
- **Pass the whole product to `weighting`.** If \\K\\ is a product of
  several size-dependent terms, build the combined product first — SSB
  uses `sweep(params@maturity, 2, params@w, "*")` (maturity \\\times\\
  body weight), yield uses fishing mortality \\\times\\ body weight.
  Bin-averaging happens inside on the weighting array as a whole, and
  the average of a product is not the product of the averages. Do not
  include `params@dw`:
  [`sizeIntegral()`](https://sizespectrum.org/mizer/reference/sizeIntegral.md)
  handles bin widths and bin-averaging automatically.
- **Extra dimensions in `weighting` are kept.** A gear × species × size
  weighting array (like
  [`getFMortGear()`](https://sizespectrum.org/mizer/reference/getFMortGear.md))
  gives a gear × species result; a weighting array whose first dimension
  is named `"time"` is lined up with the times of the simulation rather
  than multiplied out against them.

The result is already an
[`ArrayTimeBySpecies`](https://sizespectrum.org/mizer/reference/ArrayTimeBySpecies.md)
when it is one, so you inherit the whole toolkit described below —
[`plot()`](https://sizespectrum.org/mizer/reference/plot.md),
[`plot2()`](https://sizespectrum.org/mizer/reference/plot2.md),
[`plotRelative()`](https://sizespectrum.org/mizer/reference/plotRelative.md),
[`addPlot()`](https://sizespectrum.org/mizer/reference/addPlot.md) — for
free. For a quantity that keeps the size dimension, and so is not an
integral over sizes, wrap it yourself:

``` r

ArraySpeciesBySize(my_size_resolved, value_name = "My index", params = params,
                   representation = "average")
```

Use `representation = "average"` for a quantity that is a bin average
(anything integrated over a bin) and `"point"` for one sampled at the
bin boundary, such as a growth rate; the tag drives the half-bin
plotting shift.

If your indicator decomposes the encounter rate — a diet or
trophic-level style quantity — see the note on
[`encounter_kernel()`](https://sizespectrum.org/mizer/reference/encounter_kernel.md)
in the [guide to extending
mizer](https://sizespectrum.org/mizer/articles/guide-extend-mizer.md)
before pairing
[`pred_kernel()`](https://sizespectrum.org/mizer/reference/setPredKernel.md)
with
[`getEncounter()`](https://sizespectrum.org/mizer/reference/getEncounter.md).

------------------------------------------------------------------------

## Plotting mizer arrays

The arrays returned by the summary and rate functions carry a mizer
array class and have their own
[`plot()`](https://sizespectrum.org/mizer/reference/plot.md) method, so
you can visualise **any** quantity without a dedicated plot function or
custom ggplot code. They also carry a `value_name`, `type`, `units` and
their `params`, and have
[`print()`](https://sizespectrum.org/mizer/reference/print.md),
[`summary()`](https://sizespectrum.org/mizer/reference/summary.md) and
[`as.data.frame()`](https://sizespectrum.org/mizer/reference/as.data.frame.md)
methods.

| Class | Typical source | [`plot()`](https://sizespectrum.org/mizer/reference/plot.md) shows |
|----|----|----|
| [`ArrayTimeBySpecies`](https://sizespectrum.org/mizer/reference/ArrayTimeBySpecies.md) | [`getBiomass(sim)`](https://sizespectrum.org/mizer/reference/getBiomass.md), [`getSSB(sim)`](https://sizespectrum.org/mizer/reference/getSSB.md), [`getYield(sim)`](https://sizespectrum.org/mizer/reference/getYield.md), [`getN(sim)`](https://sizespectrum.org/mizer/reference/getN.md) | value vs time, one line per species |
| [`ArraySpeciesBySize`](https://sizespectrum.org/mizer/reference/ArraySpeciesBySize.md) | [`getFeedingLevel(params)`](https://sizespectrum.org/mizer/reference/getFeedingLevel.md), [`getPredMort(params)`](https://sizespectrum.org/mizer/reference/getPredMort.md), [`getEncounter(params)`](https://sizespectrum.org/mizer/reference/getEncounter.md) | value vs size, one line per species |
| [`ArrayTimeBySpeciesBySize`](https://sizespectrum.org/mizer/reference/ArrayTimeBySpeciesBySize.md) | [`getFMort(sim)`](https://sizespectrum.org/mizer/reference/getFMort.md), [`getPredMort(sim)`](https://sizespectrum.org/mizer/reference/getPredMort.md) | one time slice vs size (set with `time`) |
| [`ArrayResourceBySize`](https://sizespectrum.org/mizer/reference/ArrayResourceBySize.md) | [`NResource(params)`](https://sizespectrum.org/mizer/reference/N.md), [`finalNResource(sim)`](https://sizespectrum.org/mizer/reference/finalN.md), [`getResourceMort(params)`](https://sizespectrum.org/mizer/reference/getResourceMort.md), [`resource_rate(params)`](https://sizespectrum.org/mizer/reference/setResource.md), [`resource_capacity(params)`](https://sizespectrum.org/mizer/reference/setResource.md), [`resource_level(params)`](https://sizespectrum.org/mizer/reference/setResource.md) | resource quantity vs size |
| [`ArrayTimeByResourceBySize`](https://sizespectrum.org/mizer/reference/ArrayTimeByResourceBySize.md) | [`NResource(sim)`](https://sizespectrum.org/mizer/reference/N.md) | one time slice vs size (set with `time`) |

``` r

plot(getBiomass(sim))          # value vs time, one line per species
plot(getFeedingLevel(params))  # value vs size, one line per species
plot(getResourceMort(params))  # plankton resource mortality vs size
```

The array plots come with a small toolkit for combining and comparing
them. Every one of these has a method for every array class in the table
above:

| Function | What it does |
|----|----|
| [`addPlot()`](https://sizespectrum.org/mizer/reference/addPlot.md) | adds a compatible array as extra lines on an existing plot |
| [`plot2()`](https://sizespectrum.org/mizer/reference/plot2.md) | compares two compatible arrays (colour = species, linetype = which object) |
| [`plotRelative()`](https://sizespectrum.org/mizer/reference/plotRelative.md) | shows the relative difference `2 (y - x) / (x + y)` between two compatible arrays |
| [`plotHover()`](https://sizespectrum.org/mizer/reference/plotHover.md) | turns any of these ggplots into a hover-enabled plotly plot |

``` r

# Add another compatible array as extra lines on an existing plot
p <- plot(getBiomass(sim), species = "Cod")
addPlot(p, getBiomass(sim), species = "Herring", linetype = "dashed")

# Compare two compatible arrays
plot2(getFMort(params), getFMort(params2), "Before", "After")
plotRelative(getEGrowth(params), getEGrowth(params2))  # relative difference

plotHover(getBiomass(sim))     # interactive (hover) version of any array plot
```

[`plot2()`](https://sizespectrum.org/mizer/reference/plot2.md) and
[`plotRelative()`](https://sizespectrum.org/mizer/reference/plotRelative.md)
prepare each of their two arrays separately, using the model attached to
that array. This matters when the two models differ in the length-weight
relationship `w = a l^b`: with `size_axis = "l"` each spectrum is then
drawn at its own lengths, and a density is rescaled by its own Jacobian.
The two length grids no longer coincide when that happens, so
[`plotRelative()`](https://sizespectrum.org/mizer/reference/plotRelative.md)
interpolates both series (linearly in the logarithm of size) onto the
union of their coordinates, restricted to the range both cover. On a
weight axis, and whenever the two models agree, the grids coincide and
nothing is approximated.

The two arrays must hold the same kind of value: comparing a `"density"`
with a `"value"` is an error, because the `type` decides both the
Jacobian and the y-axis scaling and there is no pair of axes that
carries both. A differing `value_name` or differing units only warn.

### Common arguments

Most analysis and plotting functions — including
[`plot()`](https://sizespectrum.org/mizer/reference/plot.md) on an array
and the dedicated `plot…()` functions below — share these optional
arguments:

| Argument | Effect |
|----|----|
| `species` | character vector — restrict to a subset of species |
| `time_range` | numeric vector — average over this time period (plots against size) |
| `tlim` | numeric vector `c(min, max)` — restrict the time axis (plots against time) |
| `wlim`/`llim` | numeric vector `c(min, max)` — restrict the size (x) axis |
| `ylim` | numeric vector `c(min, max)` — restrict the value (y) axis |
| `highlight` | character vector — draw named species with thicker lines |
| `total` | logical — add a line for the community total. The total of *everything the object holds*, so it does not change when you select species or hide the resource; on a length axis it is summed at equal length |
| `background` | logical — whether species marked with [`markBackground()`](https://sizespectrum.org/mizer/reference/markBackground.md) are drawn. They are drawn only when the selection asks for them, always under a single grey `"Background"` legend entry; `background = FALSE` removes them |
| `log_x`, `log_y` | logical — log-scale the x or y axis |
| `size_axis` | `"w"` (default) or `"l"` — plot against weight or against length |

`wlim`/`llim` (size axis) and `ylim` (value axis) only set the **visible
window**: data outside the range is hidden but nothing is recomputed. To
change the underlying numbers — for example the size range that a
biomass is summed over — pass `min_w`/`max_w` (or `min_l`/`max_l`) to
the `get…()` function instead, e.g. `plot(getBiomass(sim, min_w = 10))`.

Which arguments apply depends on the array’s shape:

- `plot(<ArrayTimeBySpecies>)` accepts `species`, `tlim`, `total`,
  `background`, `highlight`, `log_x`, `log_y`, `ylim`.
- `plot(<ArraySpeciesBySize>)` accepts `species`, `highlight`, `total`,
  `background`, `log_x`, `log_y`, `wlim`, `llim`, `ylim`, `size_axis`,
  `per_log_size`, `all.sizes`. `size_axis`, `llim` and `per_log_size`
  belong to the size shapes only — a plot against time has no size axis
  to convert.
- `plot(<ArrayTimeBySpeciesBySize>)` takes one time slice and hands it
  to the `ArraySpeciesBySize` method, so it accepts everything that
  method does plus `time` (default: the last time step). It has no
  `tlim`: only one time is shown.
- `plot(<ArrayResourceBySize>)` accepts `log_x`, `log_y`, `wlim`,
  `llim`, `ylim`, `size_axis`, `per_log_size`. The resource is a single
  spectrum, so there is nothing for `species`, `highlight`, `total` or
  `background` to select.
- `plot(<ArrayTimeByResourceBySize>)` accepts the same as
  `ArrayResourceBySize` plus `time`.

All five also accept `return_data = TRUE`, which returns the data frame
behind the plot instead of the plot, and `y_ticks` to set the number of
y-axis ticks.

### What kind of value an array holds

Every mizer array declares what kind of quantity it holds, in its `type`
attribute, because two kinds need handling that the numbers alone do not
reveal:

| `type` | Meaning | What the plots do with it |
|----|----|----|
| `"value"` | a rate, an amount — the default | nothing special |
| `"density"` | an amount per gram of body weight | converts the values, not just the axis, when plotted against length |
| `"proportion"` | a fraction | shows the whole of the interval from 0 to 1 on a linear y axis |

Read it with
[`array_type(x)`](https://sizespectrum.org/mizer/reference/array_type.md),
and set it when you build an array of your own:

``` r

ArraySpeciesBySize(x, value_name = "Number density", units = "1/g",
                   type = "density", params = params)
```

### Plotting densities

A density is an amount *per unit size*, so its numerical value depends
on which size variable it is a density in. Changing that variable —
weight to length, or size to log size — therefore changes the plotted
**values**, not just the axis: it needs a Jacobian factor. The plot
functions apply it for you, for the arrays that declare themselves
densities:

| Source | Density |
|----|----|
| [`initialN(params)`](https://sizespectrum.org/mizer/reference/initialN-set.md), [`finalN(sim)`](https://sizespectrum.org/mizer/reference/finalN.md), [`N(sim)`](https://sizespectrum.org/mizer/reference/N.md), [`get_initial_n(params)`](https://sizespectrum.org/mizer/reference/get_initial_n.md) | consumer number density, per gram |
| [`initialNResource(params)`](https://sizespectrum.org/mizer/reference/initialNResource-set.md), [`finalNResource(sim)`](https://sizespectrum.org/mizer/reference/finalN.md), [`NResource(sim)`](https://sizespectrum.org/mizer/reference/N.md) | resource number density, per gram |
| [`resource_capacity(params)`](https://sizespectrum.org/mizer/reference/setResource.md) | resource carrying capacity, per gram |
| [`getFluxGradient(params)`](https://sizespectrum.org/mizer/reference/getFluxGradient.md) | rate of change of the flux, per gram per year |

Which density you get is set by two independent arguments: `size_axis`
chooses the size variable and `per_log_size` chooses whether the values
are per size or per logarithmic size. The factors are built from the
length-weight relationship \\w = a\\ l^b\\ of each species, taken from
the `a` and `b` columns of
[`species_params`](https://sizespectrum.org/mizer/reference/species_params.md):

| Argument                                  | Factor                   |     |
|-------------------------------------------|--------------------------|-----|
| `size_axis = "w"`, `per_log_size = TRUE`  | \\dw/d\log w = w\\       |     |
| `size_axis = "l"`, `per_log_size = FALSE` | \\dw/dl = b\\ w / l\\    |     |
| `size_axis = "l"`, `per_log_size = TRUE`  | \\dw / d\log l = b\\ w\\ |     |

**`log_x` does not change the y-axis.** Showing size on a logarithmic
axis is a display choice; you need to use `per_log_size` to convert a
density per unit size into a density per logarithmic size interval.
Conflating the two is the usual reason a spectrum looks like it has the
wrong slope.

``` r

plot(initialN(params), per_log_size = TRUE)                 # per log weight
plot(initialN(params), size_axis = "l", per_log_size = TRUE) # per log length
plot(initialNResource(params), per_log_size = TRUE)          # resource too
```

------------------------------------------------------------------------

## Plotting size spectra

[`plotSpectra()`](https://sizespectrum.org/mizer/reference/plotSpectra.md)
is the function you want for plots of the abundance or biomass density
against size, one line per species. Unlike a plain
[`plot()`](https://sizespectrum.org/mizer/reference/plot.md) of a
species density array it also overlays the resource spectrum
(`resource = TRUE`, the default). Which density it shows is set by
`biomass` and `per_log_size`, described below.

By default it shows the final time step of a simulation; pass
`time_range` to average over a period, or give it a `MizerParams` object
to see the current state. The common arguments above all apply, and
[`plotlySpectra()`](https://sizespectrum.org/mizer/reference/plotSpectra.md)
is the interactive twin.

``` r

plotSpectra(params)                                   # spectra of the current state
plotSpectra(sim, per_log_size = TRUE, time_range = 1990:2000)
plotSpectra(sim, species = c("Cod", "Herring"), resource = FALSE)
plotSpectra(sim, biomass = TRUE, size_axis = "l")     # biomass density against length
```

**The resource has its own length convention.** It is a composite of
many taxa, so instead of a taxonomic weight-length relationship it uses
the equivalent spherical diameter of an organism with the density of
water (`a = pi/6`, `b = 3`, in
[`resource_params()`](https://sizespectrum.org/mizer/reference/resource_params.md)).
It therefore appears on a length axis, but measured differently from the
fish: a fish of a given weight is about 3.7 times longer than a sphere
of that weight. That gap at the resource-consumer boundary is real
biology, not an artefact.

### Which density a spectrum plot shows

[`plotSpectra()`](https://sizespectrum.org/mizer/reference/plotSpectra.md),
[`plotSpectra2()`](https://sizespectrum.org/mizer/reference/plotSpectra2.md)
and [`animate()`](https://sizespectrum.org/mizer/reference/animate.md)
describe the plotted quantity with two independent logical arguments:

|                   | `per_log_size = FALSE` | `per_log_size = TRUE`        |
|-------------------|------------------------|------------------------------|
| `biomass = FALSE` | number density         | number density per log size  |
| `biomass = TRUE`  | biomass density        | biomass density per log size |

The older single `power` argument is the sum of the two (0, 1, 1, 2
across that table) and is still accepted.

### Cumulative distributions

[`plotCDF(object, species, biomass, normalise)`](https://sizespectrum.org/mizer/reference/plotCDF.md)
plots cumulative abundance or biomass over size — steadier than a
density spectrum for eyeballing where biomass sits. `biomass = TRUE`
(default) accumulates biomass, `biomass = FALSE` accumulates numbers;
`normalise = FALSE` plots the cumulative total rather than the
proportion. The `per_log_size` argument is not used: a cumulative total
does not depend on it.

``` r

plotCDF(NS_params, species = c("Cod", "Herring"))
plotCDF(NS_sim, biomass = FALSE, normalise = FALSE)
```

### Comparing two size distributions

| Function | Shows |
|----|----|
| [`plotSpectra2(object1, object2, name1, name2)`](https://sizespectrum.org/mizer/reference/plotSpectra2.md) | two abundance spectra overlaid |
| [`plotSpectraRelative(object1, object2)`](https://sizespectrum.org/mizer/reference/plotSpectraRelative.md) | relative difference of two spectra |
| [`plotCDF2(object1, object2, name1, name2)`](https://sizespectrum.org/mizer/reference/plotCDF2.md) | two cumulative distributions overlaid |

``` r

plotSpectra2(params, params2, "Before", "After")
plotSpectraRelative(params, params2)         # 2 (N2 - N1) / (N1 + N2)
plotCDF2(sim, sim2, "Unfished", "Fished")
```

------------------------------------------------------------------------

## Animating through time

[`animate()`](https://sizespectrum.org/mizer/reference/animate.md) plays
a spectrum or array through the course of a simulation
([`animateSpectra()`](https://sizespectrum.org/mizer/reference/animate.md)
is a retained alias).

``` r

animate(sim)                 # abundance spectra over time
animate(getFMort(sim))       # an ArrayTimeBySpeciesBySize over time
animate(NResource(sim))      # an ArrayTimeByResourceBySize over time
```

[`animate()`](https://sizespectrum.org/mizer/reference/animate.md)
accepts most of the common arguments from
[`plot()`](https://sizespectrum.org/mizer/reference/plot.md).

------------------------------------------------------------------------

## Scanning a model over a range — `scanModel()`

Everything above measures one model.
[`scanModel()`](https://sizespectrum.org/mizer/reference/scanModel.md)
measures a *family* of them: it varies one aspect of the model over a
range of values and, at each value, projects until the model settles and
measures a quantity on the attractor it settled on. The result is a
[`MizerScan`](https://sizespectrum.org/mizer/reference/MizerScan.md), a
data frame that knows how to plot itself.

``` r

scan <- scanModel(params,
                  scan_values = seq(0, 1.5, 0.1),
                  set_func    = scanFishingMortality("Cod"),  # what to vary
                  value_func  = getYield,                      # what to measure
                  species     = "Cod")
plot(scan)
attr(scan, "at_max")      # the F at which the yield is largest, i.e. F_MSY
```

- **`set_func(params, value)`** returns a modified `MizerParams`.
  Ready-made ones:
  [`scanEffort(gear)`](https://sizespectrum.org/mizer/reference/scanEffort.md)
  for fishing effort, `scanFishingMortality(species, gear)` for the
  fishing mortality on one species with the rest of the fishing left
  alone,
  [`scanSpeciesParam(species, parameter)`](https://sizespectrum.org/mizer/reference/scanEffort.md)
  for a species parameter. Any function of `(params, value)` will do, as
  long as it is **idempotent** — with `continuation = TRUE` it is handed
  the object it returned at the previous scan value, so it must set
  rather than accumulate.
- **`value_func(sim)`** returns a time-by-series object. `getBiomass`,
  `getYield`, `getSSB`, `getN` and `sizeIntegral` all work unchanged,
  and so does a plain numeric vector over time such as `getMeanWeight`.
  Give extra arguments with a closure, not through `...`:
  `value_func = function(sim) getBiomass(sim, min_w = 10)`.

Scanning something that has nothing to do with fishing needs no more
than a two-line setter, because
[`project()`](https://sizespectrum.org/mizer/reference/project.md) takes
the effort from the params object:

``` r

plot(scanModel(params,
               scan_values = 10^seq(10, 12, length.out = 9),
               set_func = function(params, value) {
                   resource_params(params)$kappa <- value
                   params
               },
               scan_name = "Resource capacity", scan_units = "g"),
     log_x = TRUE)
```

### What the band means

How the quantity is measured depends on what the model settled on, which
[`projectUntilSettled()`](https://sizespectrum.org/mizer/reference/projectUntilSettled.md)
reports:

| Attractor | What is measured |
|----|----|
| fixed point | read off the settled state, no further projection; `ymin == ymax` |
| limit cycle | averaged over **exactly one period**; `ymin`/`ymax` give the range |
| neither | averaged over `t_sample` years, and the scan values are named in a message |

So the band on the plot is the range of the oscillation, and a Hopf
bifurcation appears as the scan value at which it opens up. Averaging
over exactly one period is what keeps the curve smooth: a window that is
not a whole number of periods leaves a residue of the oscillation in the
average. If you need the average more accurately, reduce `dt` — do
**not** lengthen the window.

Points where the model settled on neither are marked with a cross, and
the `residual` column says how fast the abundances were still changing
there, in 1/year. Treat those points as provisional and raise `t_max`.

[`plot()`](https://sizespectrum.org/mizer/reference/plot.md) takes
`style = "ribbon"` (default: the average as a line inside the band),
`"envelope"` (lines along the edges, no average) or `"line"` (no band),
plus `mark_max`, `reference_lines` and the usual `log_x`/`log_y`/`ylim`
arguments. A bifurcation diagram over fishing effort is
[`scanEffort()`](https://sizespectrum.org/mizer/reference/scanEffort.md)
with `style = "envelope"`.

### `plotYieldVsF()`

The one scan common enough to have its own function.
`plotYieldVsF(params, species)` is
[`scanModel()`](https://sizespectrum.org/mizer/reference/scanModel.md)
with
[`scanFishingMortality()`](https://sizespectrum.org/mizer/reference/scanEffort.md)
and
[`getYield()`](https://sizespectrum.org/mizer/reference/getYield.md),
drawn with the peak marked, so the fishing mortality at the peak is
(F\_{MSY}):

``` r

plotYieldVsF(NS_params, "Cod", F_max = 1.5)
scan <- plotYieldVsF(NS_params, "Cod", F_max = 1.5, return_data = TRUE)
attr(scan, "at_max")      # F_MSY for Cod
```

The y axis is linear by default, because the yield is exactly zero at
`F = 0`. The current fishing mortality is drawn as a `"Current F"`
reference line. If the species already has an `F_MSY` species parameter
it is also drawn as a reference line, so the value the model gives can
be compared with the value assumed.

------------------------------------------------------------------------

## Dedicated plot functions

Besides the spectrum plots above, mizer has a dedicated `plot…()`
function for each of the common summary quantities. Each is a shortcut
for [`plot()`](https://sizespectrum.org/mizer/reference/plot.md) applied
to the matching `get…()` array
(e.g. [`plotBiomass(sim)`](https://sizespectrum.org/mizer/reference/plotBiomass.md)
is `plot(getBiomass(sim))`). They accept the common arguments above, and
each has a `plotly…()` counterpart
(e.g. [`plotlyBiomass()`](https://sizespectrum.org/mizer/reference/plotBiomass.md))
for interactive use — the array
[`plot()`](https://sizespectrum.org/mizer/reference/plot.md)s use
[`plotHover()`](https://sizespectrum.org/mizer/reference/plotHover.md)
instead.

- **Against time**: `plotBiomass(sim)`,
  [`plotYield(sim)`](https://sizespectrum.org/mizer/reference/plotYield.md)
- **Against body size** (final time step by default, or pass
  `time_range`):
  [`plotFeedingLevel(sim)`](https://sizespectrum.org/mizer/reference/plotFeedingLevel.md),
  [`plotPredMort(sim)`](https://sizespectrum.org/mizer/reference/plotPredMort.md),
  [`plotFMort(sim)`](https://sizespectrum.org/mizer/reference/plotFMort.md)
- **Distinct plots**:
  - [`plotYieldGear(sim)`](https://sizespectrum.org/mizer/reference/plotYieldGear.md)
    — yield vs time faceted by gear (one panel per gear)
  - [`plotGrowthCurves(sim)`](https://sizespectrum.org/mizer/reference/plotGrowthCurves.md)
    — size at age rather than a size spectrum
  - [`plotDiet(params)`](https://sizespectrum.org/mizer/reference/plotDiet.md)
    — stacked diet composition by prey
- **Calibration**:
  [`plotBiomassObservedVsModel(params)`](https://sizespectrum.org/mizer/reference/plotBiomassObservedVsModel.md)
  and
  [`plotYieldObservedVsModel(params)`](https://sizespectrum.org/mizer/reference/plotYieldObservedVsModel.md);
  the latter takes a `gear` argument that restricts both the modelled
  and the observed catch to the named gears. See the [guide to reaching
  steady state and
  calibrating](https://sizespectrum.org/mizer/articles/guide-calibrate-model.md).
- **Overview**: `plot(sim)` combines several panels; `plot(params)`
  shows the same panels for a model’s steady state (without the
  biomass-through-time panel).

------------------------------------------------------------------------

## Working with ggplot2

All plotting functions return a ggplot2 object, so you can customise
them:

``` r

library(ggplot2)
p <- plotBiomass(sim, species = c("Cod", "Herring"))
p + theme_bw() + labs(title = "Biomass through time")
p + geom_hline(aes(yintercept = 1e10), linetype = "dashed")
```

Species line colours and types come from the `linecolour`/`linetype`
slots of the `MizerParams`; change them there for consistent styling
across every plot:

``` r

params <- setColours(params, c("Cod" = "darkblue"))
params <- setLinetypes(params, c("Cod" = "dashed"))
```

For interactive exploration prefer the `plotly…()` twin of a named
function, or
[`plotHover()`](https://sizespectrum.org/mizer/reference/plotHover.md)
for the compositional array plots.

------------------------------------------------------------------------

## Quick reference

``` r

# ── Accessing raw arrays ───────────────────────────────────────────────────────
N(sim)                  # time × species × size
NResource(sim)          # time × size
finalN(sim)             # species × size  (last time step)
finalNResource(sim)     # size            (last time step)

# ── Species biomass / abundance / yield (time × species) ──────────────────────
getBiomass(sim)         # total biomass
getSSB(sim)             # spawning stock biomass
getN(sim)               # total abundance (numbers)
getYield(sim)           # catch in weight
getYieldGear(sim)       # catch by gear (time × gear × species)

# ── Rates at size (time × species × size) ─────────────────────────────────────
getFeedingLevel(sim)    # satiation (0 = starving, 1 = full)
getPredMort(sim)        # predation mortality
getFMort(sim)           # fishing mortality
getFMortGear(sim)       # fishing mortality by gear (time × gear × species × size)

# ── Diet and trophic (species × size × …) ────────────────────────────────────
getDiet(params)                  # proportion of diet from each prey
getTrophicLevel(params)          # trophic level at size (species × size)
getTrophicLevelBySpecies(params) # mean trophic level (species)

# ── Community indicators (time series) ────────────────────────────────────────
getProportionOfLargeFish(sim, threshold_w = 100)
getMeanWeight(sim)
getMeanLength(sim)               # needs species_params `a` and `b`
getMeanMaxWeight(sim)
getCommunitySlope(sim)          # returns data.frame with slope, intercept, R²

# ── Your own indicator: an integral over the size spectrum ────────────────────
sizeIntegral(params, weighting = w(params), min_w = 10, max_w = 5000)  # = getBiomass()
sizeIntegral(sim, weighting = sweep(maturity(params), 2, w(params), "*"),  # = getSSB()
             value_name = "SSB", units = "g")   # pass the whole product as weighting
# no dw, no bin-averaging by hand, no size-grid subsetting: sizeIntegral does it
ArraySpeciesBySize(x, params = params, representation = "average")  # size-resolved
bin_average_weight(K, params)   # the primitive, if you are not doing an integral
encounter_kernel(params)        # kernel getEncounter() uses; NOT pred_kernel()

# ── Plot any array directly, plus combine / compare tools ─────────────────────
plot(getResourceMort(params))   # any get*() array plots directly
p <- plot(getBiomass(sim), species = "Cod")
addPlot(p, getBiomass(sim), species = "Herring", linetype = "dashed")  # add lines
plot2(getFMort(params), getFMort(params2), "Before", "After")  # compare arrays
plotRelative(getEGrowth(params), getEGrowth(params2))          # relative diff
plotHover(getBiomass(sim))      # interactive (hover) version of an array plot

# ── Size spectra and other densities ──────────────────────────────────────────
plotSpectra(sim)        # abundance spectra vs size (+ resource & background)
plotCDF(sim)            # cumulative biomass/abundance over size
animate(sim)            # animate spectra through time
plotSpectra(sim, biomass = TRUE)                     # biomass rather than number
plotSpectra(sim, per_log_size = TRUE)                # density per log size
plotSpectra(sim, size_axis = "l")                    # x axis in length, not weight
plotSpectra(sim, log_x = TRUE)   # display only: does NOT change the y density

# ── Compare two simulations or models ─────────────────────────────────────────
plotSpectra2(params, params2, "Before", "After")
plotSpectraRelative(params, params2)      # relative difference of spectra
plotCDF2(sim, sim2, "Unfished", "Fished")

# ── Dedicated plot functions ──────────────────────────────────────────────────
# Each plot*() is a shortcut for plot() on the matching get*() array, and each has
# an interactive plotly*() twin (plotlyBiomass(), plotlySpectra(), …).
plot(sim)               # 5-panel summary
plotBiomass(sim)        # biomass vs time
plotYield(sim)          # yield vs time
plotYieldGear(sim)      # yield vs time, faceted by gear
plotFeedingLevel(sim)   # feeding level vs size
plotPredMort(sim)       # predation mortality vs size
plotFMort(sim)          # fishing mortality vs size
plotGrowthCurves(sim)   # size vs age
plotDiet(params, species = "Cod")  # diet composition vs size
```
