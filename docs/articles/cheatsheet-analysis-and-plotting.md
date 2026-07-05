# Cheatsheet: Analysis and Plotting

This cheatsheet gives a quick overview of the functions available in
mizer for analysing the results of simulations and creating plots. For
full documentation of each function, follow the links.

Most functions accept either a `MizerSim` object (returning a **time
series**) or a `MizerParams` object (returning a **single value** from
the initial state).

------------------------------------------------------------------------

## Accessing simulation arrays

These functions extract raw arrays from a `MizerSim` object.

| Function | Returns | Dimensions |
|----|----|----|
| [`N(sim)`](https://sizespectrum.org/mizer/reference/N.md) | species abundance density | time × species × size |
| [`NResource(sim)`](https://sizespectrum.org/mizer/reference/NResource.md) | resource abundance density | time × size |
| [`finalN(sim)`](https://sizespectrum.org/mizer/reference/finalN.md) | species abundance at last time | species × size |
| [`finalNResource(sim)`](https://sizespectrum.org/mizer/reference/finalNResource.md) | resource abundance at last time | size |
| [`getEffort(sim)`](https://sizespectrum.org/mizer/reference/getEffort.md) | fishing effort | time × gear |
| [`getTimes(sim)`](https://sizespectrum.org/mizer/reference/getTimes.md) | saved time steps | vector |

**Example:** Pull out Cod abundance:

``` r

N(sim)[, "Cod", ]          # time × size for Cod
N(sim)["2010", "Cod", ]    # size vector for Cod in year 2010
finalN(sim)["Cod", ]       # size vector for Cod at the final time step
```

------------------------------------------------------------------------

## Summary functions

These compute derived quantities from abundances. All accept `MizerSim`
or `MizerParams`. See
[`?summary_functions`](https://sizespectrum.org/mizer/reference/summary_functions.md)
for the full list.

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

**Example:**

``` r

getSSB(sim)                              # SSB of all species over time
getBiomass(sim, min_w = 10, max_w = 1e4) # biomass of 10g–10kg fish
getYield(sim)["2010", ]                  # yield in year 2010
```

The result is an `ArrayTimeBySpecies` (time × species) or
`ArraySpeciesBySize` (species × size), which can be plotted directly
with [`plot()`](https://sizespectrum.org/mizer/reference/plot.md) — see
[Plotting any array directly](#plotting-any-array-directly-with-plot)
below.

------------------------------------------------------------------------

## Indicator functions

These compute community-level indicators. All accept `MizerSim` (time
series) or `MizerParams` (single value from initial state). See
[`?indicator_functions`](https://sizespectrum.org/mizer/reference/indicator_functions.md).

| Function | Key arguments | Returns |
|----|----|----|
| [`getProportionOfLargeFish(sim)`](https://sizespectrum.org/mizer/reference/getProportionOfLargeFish.md) | `threshold_w = 100`, `biomass_proportion` | proportion of large fish through time |
| [`getMeanWeight(sim)`](https://sizespectrum.org/mizer/reference/getMeanWeight.md) | `min_w`, `max_w`, `species` | mean community weight through time |
| [`getMeanMaxWeight(sim)`](https://sizespectrum.org/mizer/reference/getMeanMaxWeight.md) | `measure = "both"/"numbers"/"biomass"` | mean asymptotic weight through time |
| [`getCommunitySlope(sim)`](https://sizespectrum.org/mizer/reference/getCommunitySlope.md) | `min_w`, `max_w`, `species` | slope, intercept, R² through time |

**Example:**

``` r

lfi <- getProportionOfLargeFish(sim, min_w = 10, max_w = 5000, threshold_w = 500)
lfi[c("1972", "2010")]

slope <- getCommunitySlope(sim, min_w = 10, max_w = 5000)
head(slope)
```

------------------------------------------------------------------------

## Plotting functions

All plotting functions return a `ggplot2` object that you can customise
further (see [Working with ggplot2](#working-with-ggplot2)). See
[`?plotting_functions`](https://sizespectrum.org/mizer/reference/plotting_functions.md).

The rest of this section starts with the general mechanism — calling
[`plot()`](https://sizespectrum.org/mizer/reference/plot.md) directly on
any array, and the arguments that control it — and only then describes
the dedicated `plot...()` functions, which are mostly shortcuts for it.

### Plotting any array directly with `plot()`

The arrays returned by the summary and rate functions carry a mizer
array class and have their own
[`plot()`](https://sizespectrum.org/mizer/reference/plot.md) method, so
you can plot any quantity without a dedicated plot function.

| Class | Typical source | [`plot()`](https://sizespectrum.org/mizer/reference/plot.md) shows |
|----|----|----|
| [`ArrayTimeBySpecies`](https://sizespectrum.org/mizer/reference/ArrayTimeBySpecies.md) | `getBiomass(sim)`, `getSSB(sim)`, `getYield(sim)`, `getN(sim)` | value vs time, one line per species |
| [`ArraySpeciesBySize`](https://sizespectrum.org/mizer/reference/ArraySpeciesBySize.md) | `getFeedingLevel(params)`, `getPredMort(params)`, `getEncounter(params)` | value vs size, one line per species |
| [`ArrayTimeBySpeciesBySize`](https://sizespectrum.org/mizer/reference/ArrayTimeBySpeciesBySize.md) | `getFMort(sim)`, `getPredMort(sim)` | one time slice vs size (set with `time`) |
| [`ArrayResourceBySize`](https://sizespectrum.org/mizer/reference/ArrayResourceBySize.md) | `NResource(params)`, `getResourceMort(params)`, `resource_rate(params)`, `resource_capacity(params)` | resource quantity vs size |

``` r

plot(getBiomass(sim))          # value vs time, one line per species
plot(getSSB(sim))
plot(getFeedingLevel(params))  # value vs size, one line per species
plot(getResourceMort(params))  # plankton resource mortality vs size
```

The array plots come with a small toolkit for combining and comparing
them:

``` r

# Add another compatible array as extra lines on an existing plot
p <- plot(getBiomass(sim), species = "Cod")
addPlot(p, getBiomass(sim), species = "Herring", linetype = "dashed")

# Compare two compatible arrays
plot2(getFMort(params), getFMort(params2), "Before", "After")
plotRelative(getEGrowth(params), getEGrowth(params2))  # relative difference

# Interactive (hover) version of any array plot
plotHover(getBiomass(sim))
```

| Function | What it does |
|----|----|
| [`addPlot()`](https://sizespectrum.org/mizer/reference/addPlot.md) | adds a compatible array as extra lines on an existing plot |
| [`plot2()`](https://sizespectrum.org/mizer/reference/plot2.md) | compares two compatible arrays (colour = species, linetype = which object) |
| [`plotRelative()`](https://sizespectrum.org/mizer/reference/plotRelative.md) | shows the relative difference between two compatible arrays |
| [`plotHover()`](https://sizespectrum.org/mizer/reference/plotHover.md) | turns any of these ggplots into a hover-enabled plotly plot |

### Common arguments

Most analysis and plotting functions — including
[`plot()`](https://sizespectrum.org/mizer/reference/plot.md) on an array
and the dedicated `plot...()` functions below — share these optional
arguments:

| Argument | Effect |
|----|----|
| `species` | character vector — restrict to a subset of species |
| `time_range` | numeric vector — average over this time period (plots against size) |
| `tlim` | numeric vector `c(min, max)` — restrict time axis (plots against time) |
| `wlim`/`llim` | numeric vector `c(min, max)` — restrict the size (x) axis (plots against size) |
| `ylim` | numeric vector `c(min, max)` — restrict the value (y) axis |
| `highlight` | character vector — draw named species with thicker lines |
| `total` | logical — add a line for the community total |
| `log_x`, `log_y`, `log` | logical — log-scale the x or y axis |

`wlim`/`llim` (size axis) and `ylim` (value axis) only set the visible
window: data outside the range is hidden but nothing is recomputed. To
change the underlying numbers — for example the size range that a
biomass is summed over — pass `min_w`/`max_w` (or `min_l`/`max_l`) to
the `get...()` function instead, e.g. `plotBiomass(sim, min_w = 10)`
restricts the calculation to fish above 10 g.

Which of these apply depends on the array’s shape:

`plot(<ArrayTimeBySpecies>)` accepts: `species`, `tlim`, `total`,
`background`, `highlight`, `log_x`, `log_y`, `ylim`.

`plot(<ArraySpeciesBySize>)` accepts: `species`, `highlight`, `log_x`,
`log_y`, `wlim`, `ylim`, `all.sizes`.

### Dedicated plot functions

Each dedicated `plot...()` function is essentially
[`plot()`](https://sizespectrum.org/mizer/reference/plot.md) applied to
the matching `get...()` array, so `plotBiomass(sim)` is
`plot(getBiomass(sim))` and `plotFeedingLevel(sim)` is
`plot(getFeedingLevel(sim))`. They accept the common arguments above.
The tables below note only where a function does something you could
**not** get by plotting the array directly.

Each also has a `plotly` counterpart
(e.g. [`plotlyBiomass()`](https://sizespectrum.org/mizer/reference/plotBiomass.md))
for interactive use — the array
[`plot()`](https://sizespectrum.org/mizer/reference/plot.md)s use
[`plotHover()`](https://sizespectrum.org/mizer/reference/plotHover.md)
instead.

#### Against time

| Function | How it relates to plotting the array directly |
|----|----|
| [`plotBiomass(sim)`](https://sizespectrum.org/mizer/reference/plotBiomass.md) | same as `plot(getBiomass(sim))` |
| [`plotYield(sim)`](https://sizespectrum.org/mizer/reference/plotYield.md) | same as `plot(getYield(sim))` |
| [`plotYieldGear(sim)`](https://sizespectrum.org/mizer/reference/plotYieldGear.md) | like [`plotYield()`](https://sizespectrum.org/mizer/reference/plotYield.md) but keeps the gear dimension, drawing one panel per fishing gear |

``` r

plotBiomass(sim, species = c("Cod", "Herring"), total = TRUE)
plotBiomass(sim, tlim = c(1980, 1990))
plotYield(sim, log_y = FALSE)
```

#### Against body size

By default these show the final time step; use `time_range` to average
over a period.

| Function | How it relates to plotting the array directly |
|----|----|
| [`plotFeedingLevel(sim)`](https://sizespectrum.org/mizer/reference/plotFeedingLevel.md) | same as `plot(getFeedingLevel(sim))` |
| [`plotPredMort(sim)`](https://sizespectrum.org/mizer/reference/plotPredMort.md) | same as `plot(getPredMort(sim))` |
| [`plotFMort(sim)`](https://sizespectrum.org/mizer/reference/plotFMort.md) | same as `plot(getFMort(sim))` |
| [`plotSpectra(sim)`](https://sizespectrum.org/mizer/reference/plotSpectra.md) | abundance/biomass spectra: additionally overlays the resource spectrum and background species, and `power` rescales the y axis |
| [`plotCDF(sim)`](https://sizespectrum.org/mizer/reference/plotCDF.md) | cumulative version of the spectrum (`normalise` for proportion vs total) |
| [`plotGrowthCurves(sim)`](https://sizespectrum.org/mizer/reference/plotGrowthCurves.md) | a distinct plot: size at age rather than a size spectrum |
| [`plotDiet(params)`](https://sizespectrum.org/mizer/reference/plotDiet.md) | a distinct plot: stacked diet composition by prey |

``` r

plotSpectra(sim, power = 2, time_range = 1990:2000)
plotFeedingLevel(sim, highlight = c("Cod", "Haddock"))
plotGrowthCurves(sim, species = "Cod", max_age = 20)
plotDiet(params, species = "Cod")
plotCDF(sim, power = 1)              # cumulative biomass; power = 0 for numbers
plotCDF(sim, normalise = FALSE)     # cumulative total rather than proportion
```

### Summary plot

``` r

plot(sim)      # 5-panel summary: feeding level, biomass, predation mort, fishing mort, spectra
plot(params)   # same panels for a model's steady state (no biomass-through-time panel)
```

### Comparing two simulations or models

These take two compatible `MizerSim` or `MizerParams` objects
(e.g. before and after a change) and show them together. For whole
spectra use the dedicated functions below; for any other rate array use
[`plot2()`](https://sizespectrum.org/mizer/reference/plot2.md) and
[`plotRelative()`](https://sizespectrum.org/mizer/reference/plotRelative.md)
from [Plotting any array
directly](#plotting-any-array-directly-with-plot).

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

### Animating spectra through time

[`animate()`](https://sizespectrum.org/mizer/reference/animate.md) plays
a spectrum or rate array through the course of a simulation.

``` r

animate(sim)                 # abundance spectra over time
animate(getFMort(sim))       # an ArrayTimeBySpeciesBySize over time
```

------------------------------------------------------------------------

### Working with ggplot2

All `plot...()` functions return a `ggplot2` object, so you can
customise them:

``` r

library(ggplot2)
p <- plotBiomass(sim, species = c("Cod", "Herring"))
p + theme_bw() + labs(title = "Biomass through time")
p + geom_hline(aes(yintercept = 1e10), linetype = "dashed")
```

Change species colours and line types via the `MizerParams` object:

``` r

params@linecolour["Cod"] <- "darkblue"
params@linetype["Cod"]   <- "dashed"
```

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
getDiet(params)                 # proportion of diet from each prey (predator × size × prey)
getTrophicLevel(params)         # trophic level at size (species × size)
getTrophicLevelBySpecies(params) # mean trophic level (species)

# ── Community indicators (time series) ────────────────────────────────────────
getProportionOfLargeFish(sim, threshold_w = 100)
getMeanWeight(sim)
getMeanMaxWeight(sim)
getCommunitySlope(sim)          # returns data.frame with slope, intercept, R²

# ── Dedicated plot functions ──────────────────────────────────────────────────
# Each plot*() is a shortcut for plot() on the matching get*() array, and each has
# an interactive plotly*() twin (plotlyBiomass(), plotlySpectra(), …).
plot(sim)               # 5-panel summary
plotBiomass(sim)        # biomass vs time
plotYield(sim)          # yield vs time
plotYieldGear(sim)      # yield vs time, faceted by gear
plotSpectra(sim)        # abundance spectra vs size (+ resource & background)
plotFeedingLevel(sim)   # feeding level vs size
plotPredMort(sim)       # predation mortality vs size
plotFMort(sim)          # fishing mortality vs size
plotGrowthCurves(sim)   # size vs age
plotDiet(params, species = "Cod")  # diet composition vs size
plotCDF(sim)            # cumulative biomass/abundance over size

# ── Plot any array directly, plus combine / compare tools ─────────────────────
plot(getResourceMort(params))   # any get*() array plots directly (resource mort vs size)
p <- plot(getBiomass(sim), species = "Cod")
addPlot(p, getBiomass(sim), species = "Herring", linetype = "dashed")  # add lines
plot2(getFMort(params), getFMort(params2), "Before", "After")  # compare arrays
plotRelative(getEGrowth(params), getEGrowth(params2))          # relative diff
plotHover(getBiomass(sim))      # interactive (hover) version of an array plot
animate(sim)                    # animate spectra through time

# ── Compare two simulations or models ─────────────────────────────────────────
plotSpectra2(params, params2, "Before", "After")
plotSpectraRelative(params, params2)      # relative difference of spectra
plotCDF2(sim, sim2, "Unfished", "Fished")
```
