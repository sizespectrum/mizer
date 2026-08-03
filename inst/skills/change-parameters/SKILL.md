---
name: change-parameters
description: >-
  Change parameters of an existing mizer model correctly. Use whenever the user
  wants to modify species parameters, size-dependent rates, fishing, the
  resource, or interactions — and especially when unsure which accessor to use:
  given_species_params() vs species_params(), changing a species parameter vs
  setting a rate array directly (setSearchVolume, setPredKernel, setParams…), or
  gear_params() vs the resource setters. Follow these rules to avoid changes
  that silently fail to propagate or get overwritten.
---

# Changing model parameters

Mizer offers several ways to change a model, and it is easy to be unsure which
one to reach for. The guiding principle:

> **Change the model at the highest level that expresses your intent, and let
> mizer propagate the change downwards.** Drop to a lower level only when you
> deliberately want to override mizer's calculation.

Every setter returns a **new** `MizerParams` — always reassign
(`params <- setResource(params, ...)`).

## The levels of a mizer model

A mizer model is built in layers, and almost every change is a choice of which
layer to reach into:

| Level | What it is | Change it with |
|---|---|---|
| **1. Size-independent parameters** | the high-level inputs: per-species scalars (`w_inf`, `beta`, `gamma`, `h`, `erepro`, …), fishing gears, resource scalars, and interactions | `species_params(params) <-`, `gear_params(params) <-`, `resource_params(params) <-`, `setInteraction()` |
| **2. Size-dependent rates** | arrays over size that mizer **calculates** from those parameters (search volume, metabolic rate, predation kernel, resource capacity, selectivity, …) | the `set…()` functions |
| **3. Rate functions** | the functions mizer calls to compute the rates during a simulation | `setRateFunction()`, `setComponent()` |

**Levels 1 and 2 cover almost all everyday work.** Most level-1 parameters exist
to *calculate* the level-2 arrays, so changing a parameter normally updates the
corresponding array automatically — that is what "propagating downwards" means.
Species parameters, gear parameters and resource parameters all sit together at
the top level: each is a set of scalar inputs from which mizer builds the arrays
below.

Level 3 is for going beyond mizer's built-in behaviour and is covered at the end.

## Level 1: species parameters

A size-structured model potentially has a huge number of parameters, because
rates must be specified at every size. Mizer assumes sensible functional forms
for the size dependence, so you only supply a small number of scalars. Mizer
also sets or calculates defaults for those you do not supply — see
[Calculation of Default Parameter Values](default_parameters.html) — and keeps
track of which values you **gave** and which it **calculated**.

| Accessor | Returns |
|---|---|
| `given_species_params(params)` | only the parameters you supplied explicitly |
| `calculated_species_params(params)` | the parameters mizer derived or defaulted |
| `species_params(params)` | everything (given, with calculated filling the gaps) |

**Rule (mizer ≥ 3.2): change species parameters with `species_params(params) <-`.**
It detects what you changed, records it as *given* (so defaults can no longer
overwrite it), and triggers recalculation of the derived scalars **and** the
size-dependent rate arrays that depend on it.

```r
species_params(params)$beta <- 150   # recorded as given; also rebuilds the predation kernel
```

`given_species_params(params) <-` does the same recording and recalculation, and
is preferable in **interactive** sessions because it additionally *warns* when
you change a parameter whose effect is overridden by another parameter you have
already given.

> **Version note.** Older guidance said to avoid `species_params(params) <-`
> because it bypassed the `given_species_params` protection and skipped
> recalculation. That was fixed in mizer 3.2. On mizer **< 3.2**, still prefer
> `given_species_params(params) <-` for edits.

Columns come back as named vectors:

```r
species_params(params)$w_mat        # named by species
given_species_params(params)$gamma  # NA where you never set it
```

For a fuller description of the individual parameters see the help page of
`species_params()`.

## How a species-parameter change propagates

Many species parameters exist only to set up a rate array; changing one re-runs
the relevant setter automatically:

| Species parameter(s) | Sets up | via |
|---|---|---|
| `gamma`, `q` | search volume | `setSearchVolume()` |
| `h`, `n` | maximum intake rate | `setMaxIntakeRate()` |
| `k`, `ks`, `p` | metabolic rate | `setMetabolicRate()` |
| `z0`, `z_ext`, `d` | external mortality | `setExtMort()` |
| `beta`, `sigma`, `pred_kernel_type` | predation kernel | `setPredKernel()` |
| `w_mat`, `w_mat25`, `w_repro_max`, `m` | reproduction allocation | `setReproduction()` |

Other species parameters are used **directly** and build no array (changing them
just changes the model): `alpha` (assimilation), `w_min` (egg size), `erepro`
and `R_max` (reproduction), `interaction_resource`, and the length–weight
parameters `a`, `b`.

## The feeding level is set by `f0`, not `h`

When `gamma` (the search-volume coefficient) is not supplied, mizer **derives it
from the target feeding level `f0`** — roughly `gamma ∝ h · f0 / (1 - f0)` — so
that the modelled steady-state feeding level equals `f0`. A consequence often
missed: raising `h` does **not** lower the feeding level, because the derived
`gamma` compensates and the feeding level stays pinned at `f0`. To make growth
more or less resource-dependent (a stronger or weaker density dependence, the
"phantom-jam" feedback), change **`f0`**, not `h`: low `f0` makes juvenile
growth strongly resource-limited, `f0` near 1 makes it nearly saturated and
resource-insensitive.

Corollary: `h = Inf` (a deliberately "no-satiation" model) makes the derived
`gamma` non-finite and throws `search_vol must not contain non-finite values` —
supply `gamma` explicitly in that case.

## Level 1: fishing gear parameters

Gears, selectivity, and catchability live in `gear_params(params)`, one row per
gear–species pair (row names `"species, gear"`). Assigning to it **does**
recompute the fishing arrays.

```r
gp <- gear_params(params)
gp["Cod, Otter", "catchability"] <- 0.8
gear_params(params) <- gp                  # rebuilds fishing mortality
```

Use `setFishing()` for supplying selectivity/catchability arrays directly or
setting baseline effort. See the `set-up-fishing` skill.

## Level 1: resource parameters

The resource works just like the species parameters.
`resource_params(params)` returns a named list of scalars (`kappa`, `lambda`,
`r_pp`, `n`, `w_pp_cutoff`) that set up the resource size-spectrum arrays, and —
like `species_params<-` — assigning to it **rebuilds those arrays** by calling
`setResource()`:

```r
resource_params(params)$kappa  <- 1e11     # rebuilds the carrying capacity (cc_pp)
resource_params(params)$lambda <- 2.05     # rebuilds cc_pp (slope)
resource_params(params)$r_pp   <- 10       # rebuilds the replenishment rate (rr_pp)
```

| Scalar(s) | Rebuilds |
|---|---|
| `kappa`, `lambda`, `w_pp_cutoff` | resource carrying capacity (`cc_pp`) |
| `r_pp`, `n` | resource replenishment rate (`rr_pp`) |

You can also set the size-resolved arrays directly. As with the species rate
arrays, doing so **freezes** them ("set manually"), so a later scalar change no
longer touches them:

| Function | Sets (and freezes) |
|---|---|
| `resource_capacity(params) <-` | the carrying capacity over size |
| `resource_rate(params) <-` | the replenishment (regeneration) rate over size |
| `resource_level(params) <-` | the resource level (fraction of capacity) |
| `setResource(params, …)` | any of the above, plus `lambda`, `n`, `w_pp_cutoff`, `resource_dynamics` |

```r
resource_capacity(params) <- my_capacity     # array over size; now frozen
resource_params(params)$kappa <- 1e11        # ignored while cc_pp is frozen
params <- setResource(params, reset = TRUE)  # unfreeze: recompute from the scalars
```

Keeping the resource at its steady state (where it replenishes at the rate at
which it is consumed) is a separate, deliberate step called *balancing*. It is
**not** triggered by `resource_params(params) <-`; it happens only in
`setResource()`, and in the array/level/dynamics setters, which balance by
default. Pass `balance = FALSE` to switch it off, e.g.
`resource_capacity(params, balance = FALSE) <- my_capacity`.

## Level 1: the interaction matrix

The species × species interaction matrix is set with `setInteraction()` and read
with `getInteraction()`. The resource interaction is the `interaction_resource`
species-parameter column instead.

```r
inter <- getInteraction(params)
inter["Cod", "Herring"] <- 0.5
params <- setInteraction(params, inter)
```

## Level 2: setting a rate array directly — and the freeze trap

Everything above changes level-1 parameters and lets mizer rebuild the level-2
arrays. This section is about reaching *into* level 2, when you need a
size-dependence mizer does not produce by default.

Each size-dependent rate has its own `set…()` function. Called with only
`params`, it **recomputes** the array from the current species parameters
(unless the array is frozen — see below):

```r
params <- setMetabolicRate(params)         # recompute from k, ks, p
params <- setParams(params)                # rebuild ALL rate arrays at once
```

Each rate array also has a **direct setter/getter** — `metab(params) <-`,
`search_vol(params) <-`, etc. (and `metab(params)`, `search_vol(params)` to read
it) — that writes an array straight in. This is equivalent to passing the array
to the `set…()` function, except the direct setter modifies `params` **in place**
while `set…()` returns a **new** object:

```r
metab(params) <- my_array                             # direct setter, in place
params <- setMetabolicRate(params, metab = my_array)  # same, returns new object
```

**Trap:** either way, supplying an array **freezes** it. Mizer marks it "set
manually" and will no longer update it from species parameters. After this,
changing the feeding species parameter has *no effect* on that rate:

```r
params <- setSearchVolume(params, search_vol = my_array)  # frozen
given_species_params(params)$gamma <- 2 * gamma           # search volume UNCHANGED now
```

To hand control back to mizer you must pass **`reset = TRUE`**. A bare
`set…(params)` will **not** recompute a frozen array — it leaves the manual
value in place and warns:

```r
params <- setSearchVolume(params, reset = TRUE)   # drop the override, recompute
```

If a user reports "I changed `gamma`/`h`/`beta` but nothing happened," suspect a
manually-set (frozen) rate array — recompute it with `set…(params, reset = TRUE)`.
The `reset` argument works the same way for every rate setter, including
`setFishing()` and `setResource()`.

The full set of size-dependent rate arrays, each with its direct setter/getter
and the `set…()` function that recomputes it:

| Rate array | Direct setter / getter | Recompute with |
|---|---|---|
| search volume | `search_vol(params) <-` | `setSearchVolume()` |
| maximum intake rate | `intake_max(params) <-` | `setMaxIntakeRate()` |
| metabolic rate | `metab(params) <-` | `setMetabolicRate()` |
| external mortality | `ext_mort(params) <-` | `setExtMort()` |
| external encounter rate | `ext_encounter(params) <-` | `setExtEncounter()` |
| external diffusion | `ext_diffusion(params) <-` | `setExtDiffusion()` |
| predation kernel | `pred_kernel(params) <-` | `setPredKernel()` |
| maturity ogive | `maturity(params) <-` | `setReproduction()` |
| reproduction allocation | `repro_prop(params) <-` | `setReproduction()` |

The fishing arrays (selectivity, catchability) and the resource arrays behave
the same way and are covered in the fishing and resource sections above.

## Level 3: replacing a rate function

Levels 1 and 2 change the *numbers* mizer uses; level 3 changes the *code*. At
every time step of a `project()` run, mizer calls a set of functions to compute
the instantaneous rates (encounter, feeding level, mortality, reproduction, …).
`setRateFunction()` swaps any one of these for a function of your own, without
touching the species parameters or the rate arrays:

```r
params <- setRateFunction(params, "Mort", "myMort")   # use myMort() for total mortality
getRateFunction(params)                               # list the current rate functions
```

The first argument names the rate to replace. The available names are the
components mizer computes internally: `Encounter`, `FeedingLevel`, `PredRate`,
`PredMort`, `FMort`, `Mort`, `EReproAndGrowth`, `ERepro`, `EGrowth`,
`Diffusion`, `ResourceMort`, `RDI`, `RDD` — and `Rates` itself, if you need to
replace the whole bundle. Each has a default `mizer…()` implementation (e.g.
`mizerMort()`) that you can call, wrap, or ignore.

Your function must accept the same arguments as the mizer function it replaces
and **return an array of the same shape**. The signature always includes the
current simulation time `t`:

```r
myRate <- function(params, n, n_pp, n_other, t, ...) { … }
```

### Time-dependent rates

This is the key reason to reach for `setRateFunction()`: **species parameters
and rate arrays are fixed for the whole simulation, but a rate *function*
receives the current time `t` and can therefore change as the simulation runs.**
That lets you express things the parameters cannot — seasonal forcing, a warming
trend, a management measure that switches on in a given year, and so on.

For example, to give total mortality an annual cycle (with `t` measured in
years), wrap the default `mizerMort()` and scale its result:

```r
seasonalMort <- function(params, t, ...) {
    mizerMort(params, t = t, ...) * (1 + 0.3 * sin(2 * pi * t))
}
params <- setRateFunction(params, "Mort", "seasonalMort")
```

or to step fishing-independent mortality up permanently from year 30 onwards:

```r
regimeShift <- function(params, t, ...) {
    factor <- if (t >= 30) 1.5 else 1
    mizerMort(params, t = t, ...) * factor
}
params <- setRateFunction(params, "Mort", "regimeShift")
```

Two practical points:

- **Extra parameters** your function needs go in `other_params(params)`, e.g.
  `other_params(params)$warming_rate <- 0.02`.
- Your functions must be defined in the **global environment or in a package** —
  mizer cannot find a function defined inside another function.

To add a whole new ecosystem component with `setComponent()`, see the
`extend-mizer` skill.

## Which should I use?

| I want to change… | Use |
|---|---|
| a per-species value (`beta`, `w_mat`, `h`, `erepro`, …) | `species_params(params) <- …` (mizer ≥ 3.2; `given_species_params(params) <-` interactively or on older mizer) |
| fishing gears / selectivity / catchability | `gear_params(params) <- …` |
| baseline effort or selectivity/catchability arrays | `setFishing(params, …)` |
| the resource (`kappa`, `lambda`, `r_pp`, …) | `resource_params(params) <- …` |
| species interactions | `setInteraction(params, …)` |
| a size-dependent rate, keeping it tied to the parameters | change the underlying species parameter |
| a size-dependent rate to a bespoke array (freezing it) | the direct setter `metab(params) <- …`, or the matching `set…()` with the array argument |
| a frozen rate back to its default form | `set…(params, reset = TRUE)` |
| everything after several edits | `setParams(params)` |
| how a rate is *computed* (e.g. to make it time-dependent) | `setRateFunction(params, …)` |
| the model's set of dynamical components | `setComponent()` (see the `extend-mizer` skill) |
