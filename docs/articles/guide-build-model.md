# Guide: Building a mizer model

This guide gives an overview of the functions used to build a mizer
model from a species parameter data frame. For bringing the model to
steady state and calibrating it, see the [guide to reaching steady state
and
calibrating](https://sizespectrum.org/mizer/articles/guide-calibrate-model.md).

Each `new…`/`set…` function **returns a new
[`MizerParams`](https://sizespectrum.org/mizer/reference/MizerParams.md)
object** — always reassign (`params <- f(params, ...)`); never modify
slots in place.

------------------------------------------------------------------------

## Choosing a constructor

| Function | Model type |
|----|----|
| [`newSingleSpeciesParams()`](https://sizespectrum.org/mizer/reference/newSingleSpeciesParams.md) | one species in a fixed background |
| [`newCommunityParams()`](https://sizespectrum.org/mizer/reference/newCommunityParams.md) | single size spectrum, no species identity |
| [`newTraitParams()`](https://sizespectrum.org/mizer/reference/newTraitParams.md) | several species differing only in asymptotic size |
| [`newMultispeciesParams()`](https://sizespectrum.org/mizer/reference/newMultispeciesParams.md) | fully general multi-species model |

Most work uses
[`newMultispeciesParams()`](https://sizespectrum.org/mizer/reference/newMultispeciesParams.md),
driven by a species parameter data frame. The rest of this material
covers that route.

------------------------------------------------------------------------

## Step 1 — Assemble the species parameters

[`species_params`](https://sizespectrum.org/mizer/reference/species_params.md)
is a data frame with **one row per species**. Only two columns are truly
required:

| Column | Meaning |
|----|----|
| `species` | the species name |
| `w_inf` | von Bertalanffy asymptotic weight (g) — the maximum-size parameter |

Mizer derives defaults for `w_max` (the computational grid boundary,
default `1.5 * w_inf`), `w_repro_max` and `w_mat` from `w_inf`.
Everything else has a sensible default or is calculated.

`w_inf` is the parameter to supply, but it need not be the one you
*have*: a table giving only `w_max`, or only lengths (`l_inf` or
`l_max`, with the length–weight parameters `a` and `b`), is accepted,
and mizer fills `w_inf` in from it and says so. Note that when `w_inf`
is taken from `w_max` the two are equal, so the `1.5 *` headroom above
the asymptotic size is not there.

Commonly supplied:

| Column | Meaning |
|----|----|
| `w_min` | Egg size (g, default 0.001) |
| `w_mat` | Maturity weight (g) |
| `beta` | Preferred predator/prey mass ratio (default 30) |
| `sigma` | Width of the lognormal predation kernel (default 2) |
| `k_vb` | von Bertalanffy K — used to derive `h` (and then `gamma`) if `h`/`gamma` absent |
| `h`, `gamma` | Max intake coefficient and search-volume coefficient (alternative to `k_vb`) |
| `a`, `b` | Length–weight conversion parameters (\\w = a l^b\\, defaults 0.01 and 3) |
| `alpha` | Assimilation efficiency (default 0.6) |
| `biomass_observed` | Observed biomass, for calibration |

Units: weights in **grams**, lengths in **cm**, time in **years**. A CSV
read with [`read.csv()`](https://rdrr.io/r/utils/read.table.html) is a
fine source; the package ships an example:

``` r

species_params <- read.csv(
    system.file("extdata", "NS_species_params.csv", package = "mizer"))
```

------------------------------------------------------------------------

## Step 2 — Create the MizerParams object

``` r

params <- newMultispeciesParams(species_params)
```

Useful optional arguments to
[`newMultispeciesParams()`](https://sizespectrum.org/mizer/reference/newMultispeciesParams.md):

| Argument | Effect |
|----|----|
| `interaction` | species × species matrix of dimensionless overlaps in `[0, 1]` (1 = full interaction, the default for every pair); scales encounter and predation mortality |
| `kappa`, `lambda`, `w_pp_cutoff` | resource spectrum coefficient, exponent, and cutoff size |
| `no_w` | number of consumer size bins on the logarithmic grid (default 100) |
| `min_w` | default egg size for species without `w_min` (default 0.001 g); the consumer grid starts at the smallest resulting species egg size |
| `max_w` | largest consumer-grid size; by default the largest species `w_max`, and it cannot be smaller than any species’ `w_max` |
| `min_w_pp` | smallest target size for the resource grid; it must be below the consumer grid and resolves to `1e-12` g when omitted |
| [`gear_params`](https://sizespectrum.org/mizer/reference/gear_params.md) | fishing gear definitions (usually omitted and configured later — see the [guide to setting up fishing](https://sizespectrum.org/mizer/articles/guide-set-up-fishing.md); defaults to a knife-edge gear catching every species) |
| [`second_order_w`](https://sizespectrum.org/mizer/reference/second_order_w.md) | use the second-order size-advection scheme; see the section “Numerical scheme: watch for numerical diffusion” in the [guide to running a mizer simulation](https://sizespectrum.org/mizer/articles/guide-run-simulation.md) |

The `no_w` bins are equally spaced in log weight between the resulting
consumer-grid minimum and `max_w`. The resource uses the same spacing
and adds smaller bins until `min_w_pp` lies in its smallest bin.

Change gears later with `gear_params(params) <- ...` or
[`setFishing()`](https://sizespectrum.org/mizer/reference/setFishing.md)
— see the [guide to setting up
fishing](https://sizespectrum.org/mizer/articles/guide-set-up-fishing.md).

------------------------------------------------------------------------

## Step 3 — Inspect what you built

``` r

summary(params)
species_params(params)     # given + calculated, one row per species
interaction_matrix(params) # the interaction matrix
gear_params(params)        # the fishing gears
resource_params(params)    # the resource scalars
```

At this point
[`newMultispeciesParams()`](https://sizespectrum.org/mizer/reference/newMultispeciesParams.md)
has given you only a **rough** initial spectrum. The model is not yet at
steady state and is not yet calibrated: that is the [guide to reaching
steady state and
calibrating](https://sizespectrum.org/mizer/articles/guide-calibrate-model.md),
which picks up from here.

------------------------------------------------------------------------

## Saving and loading a model

A finished model is worth persisting so you don’t rebuild it every
session. Use mizer’s own save/restore functions rather than bare
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html) — they store the
model in a version-stable form:

``` r

saveParams(params, "cod_model.rds")    # write a MizerParams to disk
params <- readParams("cod_model.rds")  # read it back
```

[`saveSim()`](https://sizespectrum.org/mizer/reference/saveParams.md)
and
[`readSim()`](https://sizespectrum.org/mizer/reference/saveParams.md) do
the same for a
[`MizerSim`](https://sizespectrum.org/mizer/reference/MizerSim.md)
object. If the model needs an extension package, these helpers preserve
its full S3 class while also checking and loading the packages it needs.
Bare
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html)/[`readRDS()`](https://rdrr.io/r/base/readRDS.html)
retain the class too, but skip those checks, upgrades and package
loading. See the [guide to using mizer extension
packages](https://sizespectrum.org/mizer/articles/guide-use-extension-packages.md).

Before saving, record who made the model and what it is for with
[`setMetadata()`](https://sizespectrum.org/mizer/reference/setMetadata.md).
This matters most when you share the model with others, because the
metadata travels with the object:

``` r

params <- setMetadata(params,
    title = "Celtic Sea model",
    description = "A multi-species model of the Celtic Sea fish community.",
    authors = list(list(name = "Your Name", email = "you@example.com")),
    url = "https://example.com/celtic-sea-model")
getMetadata(params)     # read the metadata back
```

All fields are optional and you can add fields of your own. mizer also
fills in `mizer_version`, `extensions`, `time_created` and
`time_modified` automatically.

------------------------------------------------------------------------

## Quick reference

``` r

# ── Build ─────────────────────────────────────────────────────────────────────
species_params <- read.csv(
    system.file("extdata", "NS_species_params.csv", package = "mizer"))
params <- newMultispeciesParams(species_params)
params <- newTraitParams()  # or newCommunityParams(), newSingleSpeciesParams()

# ── Inspect ───────────────────────────────────────────────────────────────────
summary(params)
species_params(params)      # given + calculated, one row per species
interaction_matrix(params)  # the interaction matrix
gear_params(params)         # the fishing gears
resource_params(params)     # the resource scalars

# ── Save / reload ─────────────────────────────────────────────────────────────
params <- setMetadata(params, title = "Celtic Sea model", ...)
saveParams(params, "model.rds")
params <- readParams("model.rds")
saveSim(sim, "sim.rds")
sim <- readSim("sim.rds")
```
