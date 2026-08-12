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

## Gotchas from the way defaults are calculated

Mizer fills in what you did not supply, and the derivations are chained:
`age_mat` → `h` → `gamma` and `ks`. Changing one parameter therefore rarely
changes just that one thing. The three sections below cover the traps this
creates; for the derivations themselves see
[Calculation of Default Parameter Values](default_parameters.html).

### A parameter and its derived partner cancel out

When `gamma` (the search-volume coefficient) is not supplied, mizer **derives it
from the target feeding level `f0`** — roughly `gamma ∝ h · f0 / (1 - f0)` — so
that the feeding level comes out at `f0`. A consequence often missed: raising
`h` does **not** lower the feeding level, because the derived `gamma`
compensates and the feeding level stays pinned at `f0`. To make growth more or
less resource-dependent (a stronger or weaker density dependence, the
"phantom-jam" feedback), change **`f0`**, not `h`: low `f0` makes juvenile
growth strongly resource-limited, `f0` near 1 makes it nearly saturated and
resource-insensitive.

`ks` does the same thing one level down. When it is not supplied it is derived
as `ks = fc · alpha · h · w_mat^(n - p)`, so raising `h` does **not** bring a
species closer to starvation: the metabolic coefficient scales with it and the
critical feeding level stays pinned at `fc`.

| To change… | change | not |
|---|---|---|
| the feeding level, how resource-limited growth is | `f0` | `h` |
| the critical feeding level, how close to starvation a species is | `fc` | `ks`, `h` |
| growth speed, with both feeding levels held fixed | `h` | — |

Corollary: `h = Inf` (a deliberately "no-satiation" model) makes the derived
`gamma` non-finite and throws `search_vol must not contain non-finite values` —
supply `gamma` explicitly in that case.

The chain also runs the other way, which is easy to forget. When `age_mat` is
given (or `k_vb`, from which mizer estimates it), `h` is *derived* from it:

```
h = (w_mat^(1-n) - w_min^(1-n)) / (age_mat * (1-n) * alpha * (f0 - fc))
```

So `f0` and `fc` are not local edits to the feeding level — they move `h`, and
`h` carries the change on into `gamma` and `ks`. The same goes for `w_mat`,
which is never just the maturity size: it feeds `h` through the formula above,
`ks` through the `w_mat^(n - p)` factor, and `w_mat25`. And `beta`/`sigma` feed
the `gamma` derivation, because that integrates the predation kernel.

If a user is puzzled that "changing `w_mat` changed my growth curve everywhere",
this is why. To pin a parameter against the chain, give it explicitly — a given
value is never recalculated.

### `f0` and `fc` are calibration targets, not model outputs

`gamma` is calibrated against an **idealised reference world**, not against your
model: a pure power-law resource `kappa * w^-lambda`, with all fish abundances
set to zero, and — under `defaults_edition() < 2` — with `interaction_resource`
forced to 1. So `f0` is the feeding level a species *would* have in that world.

- The **realised** feeding level in an assembled model differs, often a lot.
  In `NS_params` the feeding levels at `w_mat` run from 0.58 to 0.91 against a
  stored `f0` of 0.6. Read the realised value with `getFeedingLevel(params)`,
  and the realised critical value with `getCriticalFeedingLevel(params)`; do not
  read `f0` and `fc` out of `species_params()` and call them the model's feeding
  levels.
- Setting `interaction_resource` to anything other than 1 leaves `gamma`
  untouched under edition 1, so the reduction falls straight through to the
  realised feeding level and can starve a species outright.
- The resource scalars that enter the calibration do **not** trigger a
  re-derivation. Changing `kappa` or `lambda` with `resource_params(params) <-`
  rebuilds the resource arrays but leaves `gamma` and `q` alone, so the
  feeding level moves, and the condition `q = n + lambda - 2` that makes the
  feeding level size-independent silently breaks. After changing `lambda`,
  reset `q` and `gamma` so they are derived afresh:

```r
resource_params(params)$lambda <- 2.2
given_species_params(params)$q <- NA       # let mizer re-derive q ...
given_species_params(params)$gamma <- NA   # ... and gamma with it
```

### A default is only consulted while the column is missing

Defaults fill in `NA`s. Once a column exists in `species_params()` — which,
after the model is built, is true of every parameter — the machinery that would
have calculated it is bypassed. Two practical consequences:

**Arguments that only feed a default are dead after construction.**
`setExtMort()` computes `z0 = z0pre * w_inf^z0exp`, but only for missing `z0`
values. On a built model `z0` is always present, so `setExtMort(params, z0pre =
2)` changes nothing at all — not even with `reset = TRUE`, which clears the rate
array's comment and not the species parameter. Set the species parameter
instead:

```r
sp <- species_params(params)
z0exp <- resource_params(params)$n - 1   # the exponent setExtMort() would use
sp$z0 <- 2 * sp$w_inf^z0exp
species_params(params) <- sp
# Setting `sp$z0 <- NA` instead hands it back to the default calculation.
```

**Giving a value switches off the calculation that would have used another one.**
Mizer warns about three of these, but only from `given_species_params<-()` —
another reason to prefer it interactively:

| If you have given… | this is ignored |
|---|---|
| `gamma` | `f0` |
| `ks` | `fc` |
| `h` | `age_mat`, `k_vb` |

Two more are silent: `h` falls back to **30** when there is no growth
information at all (no `age_mat`, no `k_vb`), and with `n != p` the derivation
of `h` is only approximate — mizer says so at message level 1, which is easily
missed.

Finally, `w_inf` is not a purely dynamic parameter. Changing it re-derives
`w_mat`, `h`, `gamma` and `ks`, but the **size grid is fixed at construction**,
so enlarging `w_inf` gives repeated warnings that "the maximum weight of a
species is larger than the maximum weight of the model". To change the size
range, rebuild the model rather than edit the parameter.

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

The size-resolved resource arrays themselves can also be set directly; that is a
level-2 change and is covered under
[the resource arrays](#the-resource-arrays) below.

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

Each rate array has a **direct setter/getter** — `metab(params)` to read it and
`metab(params) <-` to write it, `search_vol(params)`/`search_vol(params) <-`,
and so on. This is the way to impose an array of your own. It must have the same
dimensions as the one it replaces (species × size, species in the model's
order), so the safest recipe is to read the current array, modify it, and write
it back:

```r
my_metab <- metab(params)                    # species × size, right dims and names
my_metab["Cod", ] <- 0.5 * w(params)^0.8     # bespoke size dependence for one species
metab(params) <- my_metab
```

Note that the direct setter modifies `params` **in place** — no reassignment.

**Trap:** supplying an array **freezes** it. Mizer marks it "set manually" and
will no longer update it from species parameters. After this, changing the
feeding species parameter has *no effect* on that rate:

```r
search_vol(params) <- my_array                    # frozen
given_species_params(params)$gamma <- 2 * gamma   # search volume UNCHANGED now
```

To hand control back to mizer, call the rate's `set…()` function with
**`reset = TRUE`**. These functions recompute the array from the current species
parameters — but a bare `set…(params)` will **not** touch a frozen array; it
leaves the manual value in place and warns:

```r
params <- setSearchVolume(params, reset = TRUE)   # drop the override, recompute
params <- setMetabolicRate(params)                # recompute, unless frozen
params <- setParams(params)                       # rebuild ALL rate arrays at once
```

Unlike the direct setters, the `set…()` functions return a **new** object, so
they must be reassigned. They also accept the array itself, as the argument
named like the direct setter:
``` r
params <- setMetabolicRate(params, metab = my_metab)
```
is equivalent to
``` r
metab(params) <- my_metab
```
and offers no advantage over it.

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

The fishing arrays (selectivity, catchability) behave the same way; see the
fishing section above and the `set-up-fishing` skill.

### The resource arrays

The resource arrays work just like the species rate arrays: you can write them
directly instead of letting `resource_params(params) <-` calculate them from the
scalars, and doing so **freezes** them ("set manually"), so a later scalar change
no longer touches them.

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

### Balancing the resource

The resource has a steady state: the abundance at which replenishment exactly
matches the rate at which consumers eat it. **Balancing** means adjusting the
replenishment rate and the carrying capacity together so that the resource
abundance the model currently holds *is* that steady state. Once a model has
been calibrated to steady state this matters a lot: an unbalanced change to the
resource leaves the resource off its steady state, and the whole model then
drifts the next time you `project()` it.

Because the balance condition ties rate and capacity together, **supply only one
of `resource_rate`, `resource_capacity` and `resource_level`** — mizer derives
the other. Supplying two of them while balancing is an error. The remaining
freedom in the pair does not change the steady state at all; it controls only
the resource *dynamics around* it:

- **high rate, low capacity** (resource level near 1): the resource regenerates
  fast and is barely depleted by consumption, so it behaves almost like a fixed
  background spectrum and competition for it has little effect;
- **low rate, high capacity** (small resource level): the resource is strongly
  depleted by consumption and the model is much more sensitive to changes in
  that competition.

`resource_level(params)` — the ratio of current abundance to capacity — is
therefore the natural single knob. Under the default semichemostat dynamics the
balance condition, in terms of that level ℓ(w), the current abundance N_R(w) and
the consumption mortality μ_R(w) from `getResourceMort()`, is just

```
r_R(w) = μ_R(w) · ℓ(w) / (1 − ℓ(w))        c_R(w) = N_R(w) / ℓ(w)
```

Balancing is a feature of `setResource()` alone. The array, level and dynamics
setters call it and so balance **by default**; `resource_params(params) <-` does
**not** balance, so changing `kappa` or `r_pp` that way rebuilds the arrays but
shifts the steady state. To change a resource coefficient and stay balanced, go
through `setResource()`:

```r
params <- setResource(params, resource_capacity = new_kappa)  # rebalances the rate
params <- setResource(params, resource_rate = new_r_pp)       # rebalances the capacity
resource_level(params) <- 0.5                                 # balances by default
resource_capacity(params, balance = FALSE) <- my_capacity     # no rebalancing
```

Two fine points:

- **Frozen arrays are protected from *incidental* balancing.** If a balancing
  operation is not given a replacement rate or capacity (e.g. you only changed
  `resource_dynamics`), a manually-set array is left alone and a warning is
  issued instead of being silently overwritten; use `reset = TRUE` to recompute
  it deliberately. But when you *do* supply one side explicitly, mizer derives
  the other even if that other one was frozen — pass `balance = FALSE` to keep it.
- **Balancing can fail**, with messages such as "I can't balance the resource if
  the capacity is less than the current abundance". That is a statement about
  the model, not a bug: no positive rate can hold the resource at an abundance
  above its capacity.

Balancing needs a `balance_<resource_dynamics>()` function; mizer supplies one
for the semichemostat and logistic dynamics. For the equations of each, and for
the full argument-by-argument behaviour, see the help pages of `setResource()`,
[`resource_semichemostat()`](../reference/resource_semichemostat.html) and
[`resource_logistic()`](../reference/resource_logistic.html).

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
