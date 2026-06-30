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
[Plotting arrays directly](#plotting-arrays-directly) below.

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

All plotting functions return a `ggplot2` object that can be further
modified. Each has a `plotly` counterpart
(e.g. [`plotlyBiomass()`](https://sizespectrum.org/mizer/reference/plotBiomass.md))
for interactive use. See
[`?plotting_functions`](https://sizespectrum.org/mizer/reference/plotting_functions.md).

### Plots against time

These show how a quantity changes over the course of the simulation.

| Function | Key arguments | Shows |
|----|----|----|
| [`plotBiomass(sim)`](https://sizespectrum.org/mizer/reference/plotBiomass.md) | `species`, `total`, `start_time`, `end_time`, `log_x`, `log_y` | total biomass per species |
| [`plotYield(sim)`](https://sizespectrum.org/mizer/reference/plotYield.md) | `species`, `total`, `log_x`, `log_y` | total yield per species |
| [`plotYieldGear(sim)`](https://sizespectrum.org/mizer/reference/plotYieldGear.md) | `species`, `total` | yield per species faceted by gear |

``` r

plotBiomass(sim)
plotBiomass(sim, species = c("Cod", "Herring"), total = TRUE)
plotBiomass(sim, start_time = 1980, end_time = 1990)
plotYield(sim, log_y = FALSE)
```

### Plots against body size

These show how a quantity varies with body size, by default at the final
time step. Use `time_range` to average over a time period.

| Function | Key arguments | Shows |
|----|----|----|
| [`plotSpectra(sim)`](https://sizespectrum.org/mizer/reference/plotSpectra.md) | `power`, `time_range`, `wlim`, `species` | abundance (or biomass) spectra |
| [`plotFeedingLevel(sim)`](https://sizespectrum.org/mizer/reference/plotFeedingLevel.md) | `time_range`, `species`, `highlight`, `log_x`, `log_y` | feeding level (0–1) at size |
| [`plotPredMort(sim)`](https://sizespectrum.org/mizer/reference/plotPredMort.md) | `time_range`, `species`, `highlight`, `log_x`, `log_y` | predation mortality at size |
| [`plotFMort(sim)`](https://sizespectrum.org/mizer/reference/plotFMort.md) | `time_range`, `species`, `highlight`, `log_x`, `log_y` | fishing mortality at size |
| [`plotGrowthCurves(sim)`](https://sizespectrum.org/mizer/reference/plotGrowthCurves.md) | `species`, `max_age`, `percentage`, `species_panel`, `log_x`, `log_y` | size at age |
| [`plotDiet(params)`](https://sizespectrum.org/mizer/reference/plotDiet.md) | `species`, `log_x`, `log_y` | diet composition at size |

``` r

plotSpectra(sim, power = 2, time_range = 1990:2000)
plotFeedingLevel(sim, highlight = c("Cod", "Haddock"))
plotGrowthCurves(sim, species = "Cod", max_age = 20)
plotDiet(params, species = "Cod")
```

### Summary plot

``` r

plot(sim)   # 5-panel summary: feeding level, biomass, predation mort, fishing mort, spectra
```

------------------------------------------------------------------------

## Plotting arrays directly

The arrays returned by summary functions carry class
`ArraySpeciesBySize` or `ArrayTimeBySpecies` and have
[`plot()`](https://sizespectrum.org/mizer/reference/plot.md) and
`ggplotly()` methods. This makes it easy to plot any quantity without a
dedicated plot function.

| Class | Typical source | [`plot()`](https://sizespectrum.org/mizer/reference/plot.md) shows |
|----|----|----|
| `ArrayTimeBySpecies` | `getBiomass(sim)`, `getSSB(sim)`, `getYield(sim)`, `getN(sim)` | value vs time, one line per species |
| `ArraySpeciesBySize` | `getFeedingLevel(params)`, `getPredMort(params)`, `getEncounter(params)` | value vs size, one line per species |

``` r

plot(getBiomass(sim))          # equivalent to plotBiomass(sim)
plot(getSSB(sim))
plot(getFeedingLevel(params))  # same as plotFeedingLevel(params)

# Add another compatible array to an existing plot
p <- plot(getBiomass(sim), species = "Cod")
addPlot(p, getBiomass(sim), species = "Herring", linetype = "dashed")

# Interactive versions
ggplotly(getBiomass(sim))
ggplotly(getEncounter(params))
```

`plot(<ArrayTimeBySpecies>)` accepts: `species`, `start_time`,
`end_time`, `total`, `background`, `highlight`, `log`, `ylim`.

`plot(<ArraySpeciesBySize>)` accepts: `species`, `highlight`, `log_x`,
`log_y`, `wlim`, `ylim`, `all.sizes`.

[`addPlot()`](https://sizespectrum.org/mizer/reference/addPlot.md) adds
an `ArrayTimeBySpecies` or `ArraySpeciesBySize` object as extra lines on
an existing compatible plot.

------------------------------------------------------------------------

## Common arguments

Most analysis and plotting functions share these optional arguments:

| Argument | Effect |
|----|----|
| `species` | character vector — restrict to a subset of species |
| `min_w`, `max_w` | restrict size range (by weight, in grams) |
| `min_l`, `max_l` | restrict size range (by length, in cm) |
| `time_range` | numeric vector — average over this time period (plot functions) |
| `start_time`, `end_time` | restrict time axis (time-series plots) |
| `highlight` | character vector — draw named species with thicker lines |
| `total` | logical — add a line for the community total |
| `log` | logical — log-scale y axis |

------------------------------------------------------------------------

## Working with ggplot2

All `plot*()` functions return a `ggplot2` object, so you can customise
them:

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

# ── Dedicated plot functions ───────────────────────────────────────────────────
plot(sim)               # 5-panel summary
plotBiomass(sim)        # biomass vs time
plotYield(sim)          # yield vs time
plotYieldGear(sim)      # yield vs time, faceted by gear
plotSpectra(sim)        # abundance spectra vs size
plotFeedingLevel(sim)   # feeding level vs size
plotPredMort(sim)       # predation mortality vs size
plotFMort(sim)          # fishing mortality vs size
plotGrowthCurves(sim)   # size vs age
plotDiet(params, species = "Cod")  # diet composition vs size

# ── Plot any summary array directly ───────────────────────────────────────────
plot(getBiomass(sim))           # ArrayTimeBySpecies → value vs time
plot(getFeedingLevel(params))   # ArraySpeciesBySize → value vs size
p <- plot(getBiomass(sim), species = "Cod")
addPlot(p, getBiomass(sim), species = "Herring", linetype = "dashed")
ggplotly(getBiomass(sim))       # interactive version

# ── Add plotly interactivity to any named plot function ───────────────────────
plotlyBiomass(sim)
plotlySpectra(sim)
plotlyFeedingLevel(sim)
# … etc.
```
