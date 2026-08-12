---
name: analyse-and-plot
description: >-
  Analyse and visualise the results of a mizer simulation or the state of a
  MizerParams object. Use whenever the user wants to extract, summarise, or plot
  size spectra, biomass, yield, SSB, abundance, feeding level, mortality, diet,
  trophic level, community indicators, growth curves, or the plankton resource —
  including comparing two simulations or animating spectra through time. Prefer
  the built-in mizer functions described here over writing custom extraction or
  ggplot code.
---

# Analysing and plotting mizer results

Mizer ships a large family of extraction, summary, and plotting functions.
**Always prefer these over hand-written array wrangling or custom ggplot code** —
they handle size-range integration, species colours/linetypes, and units for you.

Most functions accept **either** a `MizerSim` object (returning a **time
series**) **or** a `MizerParams` object (returning a **single value** from the
initial state). So `getBiomass(sim)` gives biomass over time, `getBiomass(params)`
gives biomass now.

To get the single value **at one time step** of a simulation, extract a
`MizerParams` snapshot with `finalParams(sim)` (last step), `initialParams(sim)`
(first step), or `getParams(sim, time_range = ...)` (averaged over a range) and
pass that in. Prefer this over indexing the time series with `idxFinalT(sim)`:

```r
getMeanMaxWeight(finalParams(sim))          # equilibrium value at the last step
getMeanMaxWeight(sim)[idxFinalT(sim), ]     # older equivalent
```

<!-- agent-only -->
If you need a function you don't see here, `grep` for `"plot"` or the specific
name in the bundled API index (path at the end of `MIZER-AGENTS.md`) before
writing custom code — don't read the whole file. The index gives you the name;
read the help page for the arguments.
<!-- /agent-only -->

## Accessing simulation arrays

These extract raw arrays from a `MizerSim` object.

| Function | Returns | Dimensions |
|---|---|---|
| `N(sim)` | species abundance density | time × species × size |
| `NResource(sim)` | resource abundance density | time × size |
| `finalN(sim)` | species abundance at last time | species × size |
| `finalNResource(sim)` | resource abundance at last time | size |
| `getEffort(sim)` | fishing effort | time × gear |
| `getTimes(sim)` | saved time steps | vector |

```r
N(sim)[, "Cod", ]          # time × size for Cod
N(sim)["2010", "Cod", ]    # size vector for Cod in year 2010
finalN(sim)["Cod", ]       # size vector for Cod at the final time step
```

## Summary functions

These compute derived quantities from abundances. All accept `MizerSim` or
`MizerParams`. See `?summary_functions` for the full list.

| Function | Returns | Dimensions |
|---|---|---|
| `getBiomass(sim, min_w, max_w)` | total biomass | time × species |
| `getSSB(sim)` | spawning stock biomass | time × species |
| `getN(sim, min_w, max_w)` | total abundance | time × species |
| `getYield(sim)` | total yield across gears | time × species |
| `getYieldGear(sim)` | yield by gear | time × gear × species |
| `getFeedingLevel(sim)` | feeding level at size | time × species × size |
| `getPredMort(sim)` | predation mortality at size | time × species × size |
| `getFMort(sim)` | fishing mortality at size | time × species × size |
| `getFMortGear(sim)` | fishing mortality by gear | time × gear × species × size |
| `getDiet(params)` | diet resolved by prey at size | predator × size × prey |
| `getTrophicLevel(params)` | trophic level at size | species × size |
| `getTrophicLevelBySpecies(params)` | mean trophic level per species | species |

**Size range:** `getBiomass()` and `getN()` accept `min_w`, `max_w`, `min_l`,
`max_l` to restrict the calculation to a size range.

```r
getSSB(sim)                              # SSB of all species over time
getBiomass(sim, min_w = 10, max_w = 1e4) # biomass of 10g–10kg fish
getYield(sim)["2010", ]                  # yield in year 2010
```

The result is a classed array that can be plotted directly with `plot()` — see
below.

## Indicator functions

These compute community-level indicators. All accept `MizerSim` (time series) or
`MizerParams` (single value from the initial state). See `?indicator_functions`.

| Function | Key arguments | Returns |
|---|---|---|
| `getProportionOfLargeFish(sim)` | `threshold_w = 100`, `biomass_proportion` | proportion of large fish through time |
| `getMeanWeight(sim)` | `min_w`, `max_w`, `species` | mean community weight through time |
| `getMeanMaxWeight(sim)` | `measure = "both"/"numbers"/"biomass"` | mean asymptotic weight through time |
| `getCommunitySlope(sim)` | `min_w`, `max_w`, `species` | slope, intercept, R² through time |

```r
lfi <- getProportionOfLargeFish(sim, min_w = 10, max_w = 5000, threshold_w = 500)
slope <- getCommunitySlope(sim, min_w = 10, max_w = 5000)
```

## Plotting any array directly with `plot()`

Every static plot mizer produces is a **ggplot2 object** you can extend with `+`.
The arrays returned by the summary and rate functions carry a mizer array class
and have their own `plot()` method, so you can visualise **any** quantity
without a dedicated plot function or custom ggplot code.

| Class | Typical source | `plot()` shows |
|---|---|---|
| `ArrayTimeBySpecies` | `getBiomass(sim)`, `getSSB(sim)`, `getYield(sim)`, `getN(sim)` | value vs time, one line per species |
| `ArraySpeciesBySize` | `getFeedingLevel(params)`, `getPredMort(params)`, `getEncounter(params)` | value vs size, one line per species |
| `ArrayTimeBySpeciesBySize` | `getFMort(sim)`, `getPredMort(sim)` | one time slice vs size (set with `time`) |
| `ArrayResourceBySize` | `NResource(params)`, `getResourceMort(params)`, `resource_rate(params)`, `resource_capacity(params)` | resource quantity vs size |

```r
plot(getBiomass(sim))          # value vs time, one line per species
plot(getFeedingLevel(params))  # value vs size, one line per species
plot(getResourceMort(params))  # plankton resource mortality vs size
```

The array plots come with a small toolkit for combining and comparing them:

| Function | What it does |
|---|---|
| `addPlot()` | adds a compatible array as extra lines on an existing plot |
| `plot2()` | compares two compatible arrays (colour = species, linetype = which object) |
| `plotRelative()` | shows the relative difference `2 (y - x) / (x + y)` between two compatible arrays |
| `plotHover()` | turns any of these ggplots into a hover-enabled plotly plot |

```r
# Add another compatible array as extra lines on an existing plot
p <- plot(getBiomass(sim), species = "Cod")
addPlot(p, getBiomass(sim), species = "Herring", linetype = "dashed")

# Compare two compatible arrays
plot2(getFMort(params), getFMort(params2), "Before", "After")
plotRelative(getEGrowth(params), getEGrowth(params2))  # relative difference

plotHover(getBiomass(sim))     # interactive (hover) version of any array plot
```

## Common arguments

Most analysis and plotting functions — including `plot()` on an array and the
dedicated `plot…()` functions below — share these optional arguments:

| Argument | Effect |
|---|---|
| `species` | character vector — restrict to a subset of species |
| `time_range` | numeric vector — average over this time period (plots against size) |
| `tlim` | numeric vector `c(min, max)` — restrict the time axis (plots against time) |
| `wlim`/`llim` | numeric vector `c(min, max)` — restrict the size (x) axis |
| `ylim` | numeric vector `c(min, max)` — restrict the value (y) axis |
| `highlight` | character vector — draw named species with thicker lines |
| `total` | logical — add a line for the community total |
| `log_x`, `log_y` | logical — log-scale the x or y axis |

`tlim` replaces the deprecated `start_time`/`end_time`, and `log_x`/`log_y`
replace the older single `log`.

`wlim`/`llim` (size axis) and `ylim` (value axis) only set the **visible
window**: data outside the range is hidden but nothing is recomputed. To change
the underlying numbers — for example the size range that a biomass is summed
over — pass `min_w`/`max_w` (or `min_l`/`max_l`) to the `get…()` function
instead, e.g. `plotBiomass(sim, min_w = 10)`.

Which arguments apply depends on the array's shape:

- `plot(<ArrayTimeBySpecies>)` accepts `species`, `tlim`, `total`, `background`,
  `highlight`, `log_x`, `log_y`, `ylim`.
- `plot(<ArraySpeciesBySize>)` accepts `species`, `highlight`, `log_x`, `log_y`,
  `wlim`, `ylim`, `all.sizes`.

## Dedicated plot functions

Each dedicated `plot…()` function is essentially `plot()` applied to the matching
`get…()` array, so `plotBiomass(sim)` is `plot(getBiomass(sim))`. They accept the
common arguments above, and each has a `plotly…()` counterpart (e.g.
`plotlyBiomass()`) for interactive use — the array `plot()`s use `plotHover()`
instead. See `?plotting_functions`.

**Against time:**

| Function | How it relates to plotting the array directly |
|---|---|
| `plotBiomass(sim)` | same as `plot(getBiomass(sim))` |
| `plotYield(sim)` | same as `plot(getYield(sim))` |
| `plotYieldGear(sim)` | like `plotYield()` but keeps the gear dimension, one panel per gear |

**Against body size.** By default these show the final time step; use
`time_range` to average over a period.

| Function | How it relates to plotting the array directly |
|---|---|
| `plotFeedingLevel(sim)` | same as `plot(getFeedingLevel(sim))` |
| `plotPredMort(sim)` | same as `plot(getPredMort(sim))` |
| `plotFMort(sim)` | same as `plot(getFMort(sim))` |
| `plotSpectra(sim)` | abundance/biomass spectra: additionally overlays the resource spectrum and background species, and `power` rescales the y axis |
| `plotCDF(sim)` | cumulative version of the spectrum (`normalise` for proportion vs total) |
| `plotGrowthCurves(sim)` | a distinct plot: size at age rather than a size spectrum |
| `plotDiet(params)` | a distinct plot: stacked diet composition by prey |

**Calibration:** `plotBiomassObservedVsModel(params)`,
`plotYieldObservedVsModel(params)`.

```r
plotBiomass(sim, species = c("Cod", "Herring"), total = TRUE)
plotSpectra(sim, power = 2, time_range = 1990:2000)
plotGrowthCurves(sim, species = "Cod", max_age = 20)
plotDiet(params, species = "Cod")
```

**Overview:** `plot(sim)` combines several panels; `plot(params)` shows the same
panels for a model's steady state (without the biomass-through-time panel).

## Cumulative distributions

`plotCDF(object, species, power, normalise)` plots cumulative abundance or
biomass over size — steadier than a density spectrum for eyeballing where
biomass sits. `power = 1` (default) gives biomass, `power = 0` gives numbers;
`normalise = FALSE` plots the cumulative total rather than the proportion.

```r
plotCDF(NS_params, species = c("Cod", "Herring"))
plotCDF(NS_sim, power = 0, normalise = FALSE)
```

## Comparing two simulations or models

For whole spectra use the dedicated functions below; for any other rate array
use `plot2()` and `plotRelative()` from the array toolkit above.

| Function | Shows |
|---|---|
| `plotSpectra2(object1, object2, name1, name2)` | two abundance spectra overlaid |
| `plotSpectraRelative(object1, object2)` | relative difference of two spectra |
| `plotCDF2(object1, object2, name1, name2)` | two cumulative distributions overlaid |

```r
plotSpectra2(params, params2, "Before", "After")
plotSpectraRelative(params, params2)         # 2 (N2 - N1) / (N1 + N2)
plotCDF2(sim, sim2, "Unfished", "Fished")
```

## Animating spectra through time

`animate()` plays a spectrum or rate array through the course of a simulation
(`animateSpectra()` is a retained alias).

```r
animate(sim)                 # abundance spectra over time
animate(getFMort(sim))       # an ArrayTimeBySpeciesBySize over time
```

## The plankton resource

Resource-related quantities come back as an `ArrayResourceBySize` — a numeric
vector over the size grid carrying a `value_name`, `units`, and its `params`,
with `print()`, `summary()`, `as.data.frame()`, and `plot()` methods. Producers
include `NResource(params)` / `finalNResource(sim)`, `getResourceMort(params)`,
`resource_rate(params)` (intrinsic birth rate), `resource_capacity(params)`
(carrying capacity), and `resource_level(params)`.

```r
plot(getResourceMort(params))   # resource mortality vs size
summary(NResource(params))
```

Time-resolved resource data (`NResource(sim)`) is an `ArrayTimeByResourceBySize`,
which `animate()` can play through time. To include the resource in a species
spectrum plot, pass `resource = TRUE` (supported by `plotSpectra()`, `plotCDF()`,
and friends).

## Working with ggplot2

All plotting functions return a ggplot2 object, so you can customise them:

```r
library(ggplot2)
p <- plotBiomass(sim, species = c("Cod", "Herring"))
p + theme_bw() + labs(title = "Biomass through time")
p + geom_hline(aes(yintercept = 1e10), linetype = "dashed")
```

Species line colours and types come from the `linecolour`/`linetype` slots of
the `MizerParams`; change them there for consistent styling across every plot:

```r
params <- setColours(params, c("Cod" = "darkblue"))
params <- setLinetypes(params, c("Cod" = "dashed"))
```

For interactive exploration prefer the `plotly…()` twin of a named function, or
`plotHover()` for the compositional array plots.
