# Guide: Setting up fishing

This guide gives an overview of how fishing is set up in mizer: gears,
selectivity, catchability, and effort. For full documentation of each
function, follow the links.

Fishing mortality in mizer is

\\F\_{g,i}(w) = S\_{g,i}(w)\\ Q\_{g,i}\\ E_g,\\

the product of gear \\g\\’s **selectivity** \\S\_{g,i}(w)\\ for species
\\i\\ at size \\w\\, its **catchability** \\Q\_{g,i}\\ for species
\\i\\, and its **effort** \\E_g\\. Selectivity and catchability are
configured through the
[`gear_params`](https://sizespectrum.org/mizer/reference/gear_params.md)
data frame and
[`setFishing()`](https://sizespectrum.org/mizer/reference/setFishing.md);
effort is set at run time.

------------------------------------------------------------------------

## The gear parameter data frame

One row per **gear–species combination** (a gear that catches three
species contributes three rows). Read and replace the table with
`gear_params(params)` and `gear_params(params) <- ...`. Core columns:

| Column | Meaning |
|----|----|
| `gear` | gear name; the column is required, but an `NA` entry defaults to the species name |
| `species` | species this row applies to; required |
| `sel_func` | name of the selectivity function (e.g. `"sigmoid_length"`); defaults to `"knife_edge"` |
| `catchability` | scales `F` for this gear–species pair; defaults to 1 under defaults edition 1 and 0.3 under edition 2 |

Plus **one column per parameter of the chosen `sel_func`**, named
exactly like the function’s argument (see below). Row names follow the
pattern `"species, gear"`.

There is one optional column, `yield_observed`, holding the observed
annual yield of that gear–species pair in grams per year. It is not used
by any rate calculation;
[`plotYieldObservedVsModel()`](https://sizespectrum.org/mizer/reference/plotYieldObservedVsModel.md)
sums it over the gears to compare the observations with the model, and
its `gear` argument restricts that comparison to the catch of selected
gears, for example `plotYieldObservedVsModel(params, gear = "Otter")`.

**Editing an existing gear table.** Pull it out, change what you need,
and assign it back — the assignment triggers recalculation of the
selectivity and catchability arrays:

``` r

gp <- gear_params(params)
gp["Cod, Otter", "catchability"] <- 0.8
gear_params(params) <- gp                 # assignment triggers recalculation
```

**Setting one up from scratch.** Assign a fresh data frame with one row
per gear–species combination. Both the `species` and `gear` columns must
be present, although `NA` gear entries are replaced by their species
name. Mizer fills the other core columns with the defaults shown in the
table above. You must supply the parameter columns of whatever
non-default `sel_func` you choose. Here two gears fish cod, each with
its own length-based selectivity:

``` r

gear_params(params) <- data.frame(
    gear         = c("Otter", "Beam"),
    species      = c("Cod",   "Cod"),
    sel_func     = "sigmoid_length",       # recycled to both rows
    l50          = c(25, 20),              # 50% selected at this length (cm)
    l25          = c(20, 15),              # 25% selected at this length (cm)
    catchability = 1
)
```

This replaces the whole gear table; mizer generates the
`"species, gear"` row names for you.

If each species is caught by only one gear, you may instead put the gear
columns directly in
[`species_params`](https://sizespectrum.org/mizer/reference/species_params.md)
when building the model; mizer copies them into `gear_params`. Later
edits to those `species_params` columns will **not** propagate — edit
`gear_params` after construction.

------------------------------------------------------------------------

## Selectivity functions

The selectivity function determines the size-dependence of fishing
mortality. Every selectivity function takes `w` as its first argument
and returns a value in `[0, 1]` at each size. Its other arguments must
appear as columns in `gear_params`.

| `sel_func` | Parameter column(s) | Shape |
|----|----|----|
| [`knife_edge`](https://sizespectrum.org/mizer/reference/knife_edge.md) (default) | `knife_edge_size` | step from 0 to 1 at that **weight** (default `w_mat`) |
| [`knife_edge_length`](https://sizespectrum.org/mizer/reference/knife_edge_length.md) | [`knife_edge_length`](https://sizespectrum.org/mizer/reference/knife_edge_length.md) | step from 0 to 1 at that **length** |
| [`sigmoid_length`](https://sizespectrum.org/mizer/reference/sigmoid_length.md) | `l50`, `l25` | smooth; lengths (cm) at 50% and 25% selection |
| [`double_sigmoid_length`](https://sizespectrum.org/mizer/reference/double_sigmoid_length.md) | `l50`, `l25`, `l50_right`, `l25_right` | dome-shaped (selects a length band) |
| [`sigmoid_weight`](https://sizespectrum.org/mizer/reference/sigmoid_weight.md) | `sigmoidal_weight`, `sigmoidal_sigma` | smooth transition in weight |

`sigmoid_length` is the most commonly used. You can also supply your own
function (first argument `w`, returns selectivity at size) and name it
in `sel_func`.

``` r

gp <- gear_params(params)
gp$sel_func <- "sigmoid_length"
gp$l50 <- 25            # 50% selected at 25 cm
gp$l25 <- 20            # 25% selected at 20 cm
gear_params(params) <- gp
```

------------------------------------------------------------------------

## The selectivity and catchability arrays

Behind the scenes mizer turns the `gear_params` table into two numeric
arrays, the ones that enter the fishing-mortality formula directly. You
can read them, and — when a `sel_func` cannot express the shape you need
— set them by hand.

| Function               | Returns         | Dimensions            |
|------------------------|-----------------|-----------------------|
| `catchability(params)` | \\Q\_{g,i}\\    | gear × species        |
| `selectivity(params)`  | \\S\_{g,i}(w)\\ | gear × species × size |

([`getCatchability()`](https://sizespectrum.org/mizer/reference/superseded_accessors.md)
and
[`getSelectivity()`](https://sizespectrum.org/mizer/reference/superseded_accessors.md)
are superseded aliases of these: they still work but should not be used
in new code.) Each has a matching setter that pushes an array straight
into the model (this routes through
[`setFishing()`](https://sizespectrum.org/mizer/reference/setFishing.md),
so validation still runs):

``` r

catchability(params)                       # gear × species matrix of Q
selectivity(params)["Otter", "Cod", ]      # the S curve for one gear–species pair

# Assign a custom selectivity curve that no sel_func produces
sel <- selectivity(params)
sel["Otter", "Cod", ] <- my_curve          # length = number of size bins, in [0, 1]
selectivity(params) <- sel                 # triggers recalculation via setFishing()
```

A direct assignment here is for shapes you cannot obtain from a
selectivity function. Setting an array by hand **freezes** it: mizer
marks it as manual and stops recalculating it from `gear_params`, so
later edits to the gear table leave your array untouched (you will see a
message saying so). To discard the hand-set array and rebuild from
`gear_params`, call `setFishing(params, reset = TRUE)`.

------------------------------------------------------------------------

## setFishing()

`setFishing(params, selectivity = NULL, catchability = NULL, reset = FALSE, initial_effort = NULL, ...)`
recomputes the fishing setup after a change. Assigning
`gear_params(params) <- ...` already triggers recalculation, so you
usually only call
[`setFishing()`](https://sizespectrum.org/mizer/reference/setFishing.md)
directly when supplying a `selectivity` or `catchability` **array**,
setting a baseline effort, or rebuilding from scratch:

``` r

params <- setFishing(params, initial_effort = c(Otter = 1, Beam = 0.5))
params <- setFishing(params, reset = TRUE)   # rebuild arrays from gear_params
```

------------------------------------------------------------------------

## Fishing effort

The model stores a **baseline effort** per gear, used when
[`project()`](https://sizespectrum.org/mizer/reference/project.md) is
called without an explicit `effort` argument.

| Function | Use |
|----|----|
| [`initial_effort(params)`](https://sizespectrum.org/mizer/reference/initial_effort.md) | read baseline effort (a named vector) |
| [`initial_effort(params) <-`](https://sizespectrum.org/mizer/reference/initial_effort.md) | set baseline effort |
| [`getEffort(sim)`](https://sizespectrum.org/mizer/reference/getEffort.md) | effort actually used over time in a simulation |

``` r

initial_effort(params) <- c(Industrial = 0, Pelagic = 1, Beam = 0.5, Otter = 0.5)
```

At run time,
[`project()`](https://sizespectrum.org/mizer/reference/project.md)
accepts `effort` in four forms:

``` r

project(params, effort = 1)                        # scalar: all gears, constant
project(params, effort = c(Otter = 0.5, Beam = 1)) # named vector: per gear, constant
project(params, effort = c(0.5, 1, 0, 0.5))        # vector in gear order, constant
project(params, effort = effort_array)             # time × gear array: through time
```

For a time-varying scenario, build a `time × gear` array with numeric,
increasing row names and gear column names:

``` r

gears <- names(initial_effort(params))
years <- 2010:2030
effort_array <- array(1, dim = c(length(years), length(gears)),
                      dimnames = list(time = years, gear = gears))
effort_array[as.character(2020:2030), "Otter"] <- 1.5   # ramp one gear from 2020
sim <- project(params, effort = effort_array)
```

------------------------------------------------------------------------

## Catchability sets the units of effort

Because \\F = S\\Q\\E\\, only the product \\Q\\E\\ is pinned down by the
fishing mortality — so **catchability defines what one unit of effort
means**. This gives two common conventions:

- **Effort = fishing mortality rate.** Set `catchability = 1`; then an
  effort of \\E\\ produces \\F = E\\ on fully-selected sizes. Simplest
  when you want to drive the model directly with F values, for example
  from a stock assessment.
- **Effort in real-world units** (vessel-days, kW-days, …). Fold the
  conversion into catchability: \\Q\\ is the fraction of the selected
  stock taken per unit of whatever effort you supply. A useful trick is
  to set each species’ catchability to its \\F\\ in a chosen reference
  year, so an effort of 1 reproduces that year’s fishing mortality and
  other years’ efforts are relative to it.

Either way, if you rescale catchability you must rescale effort
inversely to keep the same \\F\\. This is why comparing modelled to
observed yields with
[`plotYieldObservedVsModel()`](https://sizespectrum.org/mizer/reference/plotYieldObservedVsModel.md)
depends on the fishing setup being fixed first.

------------------------------------------------------------------------

## Inspecting the fishing setup

``` r

gear_params(params)        # the gear table
catchability(params)       # Q array (gear × species)
selectivity(params)        # S array (gear × species × size)
initial_effort(params)     # baseline effort per gear
plotFMort(params)          # realised fishing mortality at size
getFMort(params)           # F by species × size
getFMortGear(params)       # F by gear × species × size
getYieldGear(sim)          # yield by gear (from a MizerSim)
```

[`getFMort()`](https://sizespectrum.org/mizer/reference/getFMort.md) and
[`plotFMort()`](https://sizespectrum.org/mizer/reference/plotFMort.md)
**sum over gears**, so a gear with a tiny catchability — a survey gear
you fit against but do not want to fish with — is invisible in them.
Inspect that gear with
[`getFMortGear()`](https://sizespectrum.org/mizer/reference/getFMortGear.md),
slicing out its own mortality:

``` r

getFMortGear(params)["Survey", , ]   # species × size, for one gear only
```

------------------------------------------------------------------------

## Quick reference

``` r

# ── Gears and selectivity ─────────────────────────────────────────────────────
gp <- gear_params(params)
gp$sel_func <- "sigmoid_length"
gp$l50 <- 25; gp$l25 <- 20
gp$catchability <- 1
gear_params(params) <- gp

# ── Effort ────────────────────────────────────────────────────────────────────
initial_effort(params) <- c(Otter = 0.5, Beam = 1)   # baseline
sim <- project(params, effort = 1)                   # constant during run
sim <- project(params, effort = effort_array)        # time × gear array

# ── Inspect ───────────────────────────────────────────────────────────────────
initial_effort(params)
plotFMort(params)
getFMortGear(params)
```
