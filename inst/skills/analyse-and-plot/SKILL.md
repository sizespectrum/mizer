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
| `getTimes(sim)` | saved time steps | time |

```r
N(sim)[, , 1]              # time × species in smallest size class
N(sim)["2010", "Cod", ]    # size vector for Cod in year 2010
finalN(sim)["Cod", ]       # size vector for Cod at the final time step
```

## Summary functions

These functions compute derived quantities from abundances. All accept `MizerSim` or
`MizerParams`. 
The result is a classed array that can be plotted directly with `plot()` — see
below.

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

## Writing your own indicator

First check that a built-in does not already cover it: most custom indicators
turn out to be `getBiomass()`/`getN()` over a size range, or one of the four
above with different arguments. If none fits, an indicator is an integral over
the size spectrum, $\int N_i(w)\, K_i(w)\, dw$, and **`sizeIntegral()` does that
integral for you**:

```r
# Biomass of each species between 10g and 5kg -- what getBiomass() does
sizeIntegral(params, weight = params@w, min_w = 10, max_w = 5000)

# The same through time, wrapped ready to plot
sizeIntegral(sim, weight = params@w, value_name = "Biomass", units = "g")
```

Give it the abundance and the weight $K$; it selects the size range, uses the
quadrature scheme the model is actually on and wraps the result in the
appropriate mizer array class. Doing the sum by hand instead means getting all
of that right yourself, silently and only for some users. Three things are worth
knowing:

- **The size range is an argument, not a subsetting step.** `min_w`/`max_w` or
  `min_l`/`max_l` are passed through to `get_size_range_array()`, which does the
  length-weight conversion per species and accepts either a single number or one
  value per species. Never subset the size grid by hand.
- **Pass the whole product as the weight.** If $K$ is a product of several
  size-dependent terms, build the product first — SSB uses
  `sweep(params@maturity, 2, params@w, "*")` as one weight, yield uses `F * w`.
  Bin-averaging happens inside, on the weight as a whole, and the average of a
  product is not the product of the averages. Do not include `params@dw` and do
  not call `bin_average_weight()` yourself: `sizeIntegral()` does both.
- **Extra dimensions of the weight are kept.** A gear × species × size weight
  (like `getFMortGear()`) gives a gear × species result; a weight whose first
  dimension is named `"time"` is lined up with the times of the simulation
  rather than multiplied out against them.

The result is already an `ArrayTimeBySpecies` when it is one, so you inherit the
whole toolkit above — `plot()`, `plot2()`, `plotRelative()`, `addPlot()` — for
free. For a quantity that keeps the size dimension, and so is not an integral
over sizes, wrap it yourself:

```r
ArraySpeciesBySize(my_size_resolved, value_name = "My index", params = params,
                   representation = "average")
```

Use `representation = "average"` for a quantity that is a bin average (anything
integrated over a bin) and `"point"` for one sampled at the bin boundary, such as
a growth rate; the tag drives the half-bin plotting shift.

If you do need the weight itself for something other than an integral over
sizes, `bin_average_weight(K, params)` is the underlying primitive: it applies
the model's quadrature scheme to `K` alone and is gated on `second_order_w()`,
so it leaves `K` untouched on the default scheme. Never bin-average the
abundance `N` or the bin widths `dw` — `N` is already a bin average and `dw` is
exact, so averaging either counts the same integral twice.

If your indicator decomposes the encounter rate — a diet or trophic-level style
quantity — see the note on `encounter_kernel()` in the `extend-mizer` skill
before pairing `pred_kernel()` with `getEncounter()`.

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

The array plots come with a small toolkit for combining and comparing them.
Every one of these has a method for every array class in the table above:

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
| `size_axis` | `"w"` (default) or `"l"` — plot against weight or against length |

`tlim` replaces the deprecated `start_time`/`end_time`, and `log_x`/`log_y`
replace the older single `log`.

`wlim`/`llim` (size axis) and `ylim` (value axis) only set the **visible
window**: data outside the range is hidden but nothing is recomputed. To change
the underlying numbers — for example the size range that a biomass is summed
over — pass `min_w`/`max_w` (or `min_l`/`max_l`) to the `get…()` function
instead, e.g. `plot(getBiomass(sim, min_w = 10))`.

`size_axis = "l"` converts the axis with the length–weight parameters `a` and
`b`, so it is unavailable for the resource, which has no species to take them
from. For a *density* it converts the y-axis too, via the appropriate Jacobian —
see the next section.

### Which density a spectrum plot shows

`plotSpectra()`, `plotSpectra2()`, `plotCDF()`, `plotCDF2()` and `animate()`
describe the plotted quantity with two independent logical arguments, each of
which contributes one factor of the weight:

| | `per_log_size = FALSE` | `per_log_size = TRUE` |
|---|---|---|
| `biomass = FALSE` | number density | number density per log size |
| `biomass = TRUE` | biomass density | biomass density per log size |

The older single `power` argument is the sum of the two (0, 1, 1, 2 across that
table) and is still accepted, but it cannot tell the two `power = 1` cells
apart — it is read as the biomass density with respect to weight, which is what
picks the y-axis label and the length-axis Jacobian. Supplying `power` together
with a flag that contradicts it is an error, so express the choice with the
flags. `plotCDF()` accepts only `per_log_size = FALSE`: a cumulative total does
not depend on the density it was accumulated from.

**`log_x` does not change the y-axis.** Showing weight on a logarithmic axis is
a display choice; converting a density per unit weight into a density per
logarithmic weight interval is `per_log_size`. Conflating the two is the usual
reason a spectrum looks like it has the wrong slope.

Which arguments apply depends on the array's shape:

- `plot(<ArrayTimeBySpecies>)` accepts `species`, `tlim`, `total`, `background`,
  `highlight`, `log_x`, `log_y`, `ylim`.
- `plot(<ArraySpeciesBySize>)` accepts `species`, `highlight`, `log_x`, `log_y`,
  `wlim`, `llim`, `ylim`, `size_axis`, `all.sizes`. `size_axis` and `llim`
  belong to this shape only — a plot against time has no size axis to convert.

## Dedicated plot functions

Each dedicated `plot…()` function is essentially `plot()` applied to the matching
`get…()` array, so `plotBiomass(sim)` is `plot(getBiomass(sim))`. They accept the
common arguments above, and each has a `plotly…()` counterpart (e.g.
`plotlyBiomass()`) for interactive use — the array `plot()`s use `plotHover()`
instead.

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
| `plotSpectra(sim)` | abundance/biomass spectra: additionally overlays the resource spectrum and background species, and `biomass`/`per_log_size` choose the plotted density (see [above](#which-density-a-spectrum-plot-shows)) |
| `plotCDF(sim)` | cumulative version of the spectrum (`normalise` for proportion vs total) |
| `plotGrowthCurves(sim)` | a distinct plot: size at age rather than a size spectrum |
| `plotDiet(params)` | a distinct plot: stacked diet composition by prey |

**Calibration:** `plotBiomassObservedVsModel(params)` and
`plotYieldObservedVsModel(params)`; the latter takes a `gear` argument that
restricts both the modelled and the observed catch to the named gears. See the
`calibrate-model` skill.

```r
plotBiomass(sim, species = c("Cod", "Herring"), total = TRUE)
plotSpectra(sim, per_log_size = TRUE, time_range = 1990:2000)
plotGrowthCurves(sim, species = "Cod", max_age = 20)
plotDiet(params, species = "Cod")
```

**Overview:** `plot(sim)` combines several panels; `plot(params)` shows the same
panels for a model's steady state (without the biomass-through-time panel).

## Cumulative distributions

`plotCDF(object, species, biomass, normalise)` plots cumulative abundance or
biomass over size — steadier than a density spectrum for eyeballing where
biomass sits. `biomass = TRUE` (default) accumulates biomass, `biomass = FALSE`
accumulates numbers; `normalise = FALSE` plots the cumulative total rather than
the proportion. Unlike in `plotSpectra()`, only `per_log_size = FALSE` is
accepted: the integral does not depend on it.

```r
plotCDF(NS_params, species = c("Cod", "Herring"))
plotCDF(NS_sim, biomass = FALSE, normalise = FALSE)
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
animate(NResource(sim))      # an ArrayTimeByResourceBySize over time
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

The array toolkit works on the resource too:

```r
plot2(resource_capacity(params), resource_capacity(params2), "Before", "After")
plotRelative(resource_capacity(params), resource_capacity(params2))
addPlot(plot(NResource(params)), resource_capacity(params))
```

A resource array holds a single spectrum, so `species`, `total` and `background`
do nothing there and warn if you set them, and `size_axis = "l"` is unavailable
because the weight-length relationship is a species parameter.

Time-resolved resource data (`NResource(sim)`) is an `ArrayTimeByResourceBySize`,
which `animate()` can play through time, and which the comparison functions above
slice at a chosen `time`. To include the resource in a species spectrum plot,
pass `resource = TRUE` (supported by `plotSpectra()`, `plotCDF()`, and friends).

```r
animate(NResource(sim))                       # resource spectrum over time
plot2(NResource(sim), NResource(sim2), time = 1990)
```

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
