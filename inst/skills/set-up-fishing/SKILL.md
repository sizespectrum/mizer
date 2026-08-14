---
name: set-up-fishing
description: >-
  Set up or change fishing in a mizer model — gears, selectivity curves,
  catchability, and effort. Use whenever the user wants to define fishing gears,
  choose or configure a selectivity function (knife_edge, sigmoid_length,
  double_sigmoid_length, sigmoid_weight), set which gear catches which species,
  change catchability, or set the fishing effort with setFishing() and the
  gear_params data frame.
---

# Setting up fishing

Fishing mortality in mizer is

$$F_{g,i}(w) = S_{g,i}(w)\, Q_{g,i}\, E_g,$$

the product of gear $g$'s **selectivity** $S_{g,i}(w)$ for species $i$ at size
$w$, its **catchability** $Q_{g,i}$ for species $i$, and its **effort** $E_g$.
Selectivity and catchability are configured through the `gear_params` data
frame and `setFishing()`; effort is set at run time.

## The gear parameter data frame

One row per **gear–species combination** (a gear that catches three species
contributes three rows). Read and replace the table with `gear_params(params)`
and `gear_params(params) <- ...`. Required columns:

| Column | Meaning |
|---|---|
| `gear` | gear name |
| `species` | species this row applies to |
| `sel_func` | name of the selectivity function (e.g. `"sigmoid_length"`) |
| `catchability` | scales `F` for this gear–species pair (default 1) |

Plus **one column per parameter of the chosen `sel_func`**, named exactly like
the function's argument (see below). Row names follow the pattern
`"species, gear"`.

**Editing an existing gear table.** Pull it out, change what you need, and
assign it back — the assignment triggers recalculation of the selectivity and
catchability arrays:

```r
gp <- gear_params(params)
gp["Cod, Otter", "catchability"] <- 0.8
gear_params(params) <- gp                 # assignment triggers recalculation
```

**Setting one up from scratch.** Assign a fresh data frame with one row per
gear–species combination. Only `species` is strictly required: `gear` defaults
to the species name, `sel_func` to `knife_edge`, `catchability` to 1, and the
`knife_edge` cut-off to `w_mat`. You must, however, supply the parameter columns
of whatever `sel_func` you choose. Here two gears fish cod, each with its own
length-based selectivity:

```r
gear_params(params) <- data.frame(
    gear         = c("Otter", "Beam"),
    species      = c("Cod",   "Cod"),
    sel_func     = "sigmoid_length",       # recycled to both rows
    l50          = c(25, 20),              # 50% selected at this length (cm)
    l25          = c(20, 15),              # 25% selected at this length (cm)
    catchability = 1
)
```

This replaces the whole gear table; mizer generates the `"species, gear"` row
names for you.

If each species is caught by only one gear, you may instead put the gear columns
directly in `species_params` when building the model; mizer copies them into
`gear_params`. Later edits to those `species_params` columns will **not**
propagate — edit `gear_params` after construction.

## Selectivity functions

The selectivity function determines the size-dependence of fishing mortality.
Every selectivity function takes `w` as its first argument and returns a value
in `[0, 1]` at each size. Its other arguments must appear as columns in
`gear_params`.

| `sel_func` | Parameter column(s) | Shape |
|---|---|---|
| `knife_edge` (default) | `knife_edge_size` | step from 0 to 1 at that **weight** (default `w_mat`) |
| `knife_edge_length` | `knife_edge_length` | step from 0 to 1 at that **length** |
| `sigmoid_length` | `l50`, `l25` | smooth; lengths (cm) at 50% and 25% selection |
| `double_sigmoid_length` | `l50`, `l25`, `l50_right`, `l25_right` | dome-shaped (selects a length band) |
| `sigmoid_weight` | `sigmoidal_weight`, `sigmoidal_sigma` | smooth transition in weight |

`sigmoid_length` is the most commonly used. You can also supply your own
function (first argument `w`, returns selectivity at size) and name it in
`sel_func`.

```r
gp <- gear_params(params)
gp$sel_func <- "sigmoid_length"
gp$l50 <- 25            # 50% selected at 25 cm
gp$l25 <- 20            # 25% selected at 20 cm
gear_params(params) <- gp
```

## The selectivity and catchability arrays

Behind the scenes mizer turns the `gear_params` table into two numeric arrays,
the ones that enter the fishing-mortality formula directly. You can read them,
and — when a `sel_func` cannot express the shape you need — set them by hand.

| Function | Returns | Dimensions |
|---|---|---|
| `catchability(params)` | $Q_{g,i}$ | gear × species |
| `selectivity(params)` | $S_{g,i}(w)$ | gear × species × size |

(`getCatchability()` and `getSelectivity()` are deprecated aliases of these.)
Each has a matching setter that pushes an array straight into the model (this
routes through `setFishing()`, so validation still runs):

```r
catchability(params)                       # gear × species matrix of Q
selectivity(params)["Otter", "Cod", ]      # the S curve for one gear–species pair

# Assign a custom selectivity curve that no sel_func produces
sel <- selectivity(params)
sel["Otter", "Cod", ] <- my_curve          # length = number of size bins, in [0, 1]
selectivity(params) <- sel                 # triggers recalculation via setFishing()
```

A direct assignment here is for shapes you cannot obtain from a selectivity
function. Setting an array by hand **freezes** it: mizer marks it as manual and
stops recalculating it from `gear_params`, so later edits to the gear table
leave your array untouched (you will see a message saying so). To discard the
hand-set array and rebuild from `gear_params`, call
`setFishing(params, reset = TRUE)`.

## setFishing()

`setFishing(params, selectivity = NULL, catchability = NULL, reset = FALSE,
initial_effort = NULL, ...)` recomputes the fishing setup after a change.
Assigning `gear_params(params) <- ...` already triggers recalculation, so you
usually only call `setFishing()` directly when supplying a `selectivity` or
`catchability` **array**, setting a baseline effort, or rebuilding from scratch:

```r
params <- setFishing(params, initial_effort = c(Otter = 1, Beam = 0.5))
params <- setFishing(params, reset = TRUE)   # rebuild arrays from gear_params
```

## Fishing effort

The model stores a **baseline effort** per gear, used when `project()` is called
without an explicit `effort` argument.

| Function | Use |
|---|---|
| `initial_effort(params)` | read baseline effort (a named vector) |
| `initial_effort(params) <-` | set baseline effort |
| `getEffort(sim)` | effort actually used over time in a simulation |

```r
initial_effort(params) <- c(Industrial = 0, Pelagic = 1, Beam = 0.5, Otter = 0.5)
```

At run time, `project()` accepts `effort` in four forms:

```r
project(params, effort = 1)                        # scalar: all gears, constant
project(params, effort = c(Otter = 0.5, Beam = 1)) # named vector: per gear, constant
project(params, effort = c(0.5, 1, 0, 0.5))        # vector in gear order, constant
project(params, effort = effort_array)             # time × gear array: through time
```

For a time-varying scenario, build a `time × gear` array with numeric,
increasing row names and gear column names:

```r
gears <- names(initial_effort(params))
years <- 2010:2030
effort_array <- array(1, dim = c(length(years), length(gears)),
                      dimnames = list(time = years, gear = gears))
effort_array[as.character(2020:2030), "Otter"] <- 1.5   # ramp one gear from 2020
sim <- project(params, effort = effort_array)
```

## Catchability sets the units of effort

Because $F = S\,Q\,E$, only the product $Q\,E$ is pinned down by the fishing
mortality — so **catchability defines what one unit of effort means**. This
gives two common conventions:

- **Effort = fishing mortality rate.** Set `catchability = 1`; then an effort of
  $E$ produces $F = E$ on fully-selected sizes. Simplest when you want to drive
  the model directly with F values, for example from a stock assessment.
- **Effort in real-world units** (vessel-days, kW-days, …). Fold the conversion
  into catchability: $Q$ is the fraction of the selected stock taken per unit of
  whatever effort you supply. A useful trick is to set each species'
  catchability to its $F$ in a chosen reference year, so an effort of 1
  reproduces that year's fishing mortality and other years' efforts are relative
  to it.

Either way, if you rescale catchability you must rescale effort inversely to
keep the same $F$. This is why yield calibration with `calibrateYield()` depends
on the fishing setup being fixed first.

## Inspecting the fishing setup

```r
gear_params(params)        # the gear table
catchability(params)       # Q array (gear × species)
selectivity(params)        # S array (gear × species × size)
initial_effort(params)     # baseline effort per gear
plotFMort(params)          # realised fishing mortality at size
getFMort(params)           # F by species × size
getFMortGear(params)       # F by gear × species × size
getYieldGear(sim)          # yield by gear (from a MizerSim)
```

`getFMort()` and `plotFMort()` **sum over gears**, so a gear with a tiny
catchability — a survey gear you fit against but do not want to fish with — is
invisible in them. Inspect that gear with `getFMortGear()`, slicing out its own
mortality:

```r
getFMortGear(params)["Survey", , ]   # species × size, for one gear only
```
