---
name: build-multispecies-model
description: >-
  Build a multi-species mizer model from a species-parameter data frame. Use
  whenever the user wants to create a MizerParams object with
  newMultispeciesParams() (or newTraitParams, newCommunityParams,
  newSingleSpeciesParams), set up an interaction matrix or fishing gears, choose
  the size grid, or save and reload a finished model. Follow this ordered
  workflow rather than guessing at parameters or writing the dynamics by hand.
  Once the object exists, bringing it to steady state and calibrating it to data
  is covered by the calibrate-model skill.
---

# Building a multi-species mizer model

Each `new…`/`set…` function **returns a new
`MizerParams` object** — always reassign (`params <- f(params, ...)`); never
modify slots in place.

## Choosing a constructor

| Function | Model type |
|---|---|
| `newSingleSpeciesParams()` | one species in a fixed background |
| `newCommunityParams()` | single size spectrum, no species identity |
| `newTraitParams()` | several species differing only in asymptotic size |
| `newMultispeciesParams()` | fully general multi-species model |

Most work uses `newMultispeciesParams()`, driven by a species parameter data
frame. The rest of this material covers that route.

## Step 1 — Assemble the species parameters

`species_params` is a data frame with **one row per species**. Only two columns
are truly required:

| Column | Meaning |
|---|---|
| `species` | the species name |
| `w_inf` | von Bertalanffy asymptotic weight (g) — the required maximum-size parameter |

Mizer derives defaults for `w_max` (the computational grid boundary, default `1.5 * w_inf`),
`w_repro_max` and `w_mat` from `w_inf`. Everything else has a
sensible default or is calculated. Commonly supplied:

| Column | Meaning |
|---|---|
| `w_mat` | Maturity weight (g) |
| `beta` | Preferred predator/prey mass ratio (default ~100) |
| `sigma` | Width of the lognormal predation kernel (default ~1.3) |
| `k_vb` | von Bertalanffy K — used to derive `h` (and then `gamma`) if `h`/`gamma` absent |
| `h`, `gamma` | Max intake coefficient and search-volume coefficient (alternative to `k_vb`) |
| `alpha` | Assimilation efficiency (default 0.6) |
| `biomass_observed` | Observed biomass, for calibration |

Units: weights in **grams**, lengths in **cm**, time in **years**. A CSV read
with `read.csv()` is a fine source; the package ships an example:

```r
species_params <- read.csv(
    system.file("extdata", "NS_species_params.csv", package = "mizer"))
```

## Step 2 — Create the MizerParams object

```r
params <- newMultispeciesParams(species_params)
```

Useful optional arguments to `newMultispeciesParams()`:

| Argument | Effect |
|---|---|
| `interaction` | species × species matrix of dimensionless overlaps in `[0, 1]` (1 = full interaction, the default for every pair); scales encounter and predation mortality |
| `gear_params` | fishing gear definitions — columns `gear`, `species`, `sel_func`, `catchability` and the selectivity parameters. Omit it and mizer builds a default knife-edge gear catching every species |
| `no_w`, `min_w`, `max_w` | size-grid resolution and range (`no_w = 100` by default) |
| `second_order_w` | use the second-order size-advection scheme; see the section "Numerical scheme: watch for numerical diffusion" in the `run-simulation` skill |

Change gears later with `gear_params(params) <- ...` or `setFishing()` — see the
`set-up-fishing` skill.

## Step 3 — Inspect what you built

```r
summary(params)
species_params(params)     # given + calculated, one row per species
getInteraction(params)     # the interaction matrix
gear_params(params)        # the fishing gears
resource_params(params)    # the resource scalars
```

At this point `newMultispeciesParams()` has given you only a **rough** initial
spectrum. The model is not yet at steady state and is not yet calibrated: that
is the `calibrate-model` skill, which picks up from here.

## Saving and loading a model

A finished model is worth persisting so you don't rebuild it every session. Use
mizer's own save/restore functions rather than bare `saveRDS()` — they store the
model in a version-stable form:

```r
saveParams(params, "cod_model.rds")    # write a MizerParams to disk
params <- readParams("cod_model.rds")  # read it back
```

`saveSim()` and `readSim()` do the same for a `MizerSim` object.

Before saving, record who made the model and what it is for with
`setMetadata()`. This matters most when you share the model with others, because
the metadata travels with the object:

```r
params <- setMetadata(params,
    title = "Celtic Sea model",
    description = "A multi-species model of the Celtic Sea fish community.",
    authors = list(list(name = "Your Name", email = "you@example.com")),
    url = "https://example.com/celtic-sea-model")
getMetadata(params)     # read the metadata back
```

All fields are optional and you can add fields of your own. mizer also fills in
`mizer_version`, `extensions`, `time_created` and `time_modified` automatically.

## Common pitfalls

- Forgetting to reassign the return value (`steady(params)` without `params <-`)
  silently discards the work.
- Skipping `steady()` after a `match…`/`calibrate…` step leaves the model off
  its steady state.
- Editing a species parameter on a bare data frame instead of through
  `species_params(params) <-` means no dependent quantity is recalculated. See
  the `change-parameters` skill.
- Passing `max_w = w_inf` to `newMultispeciesParams()` errors, because `w_max`
  still defaults to `1.5 * w_inf` and the grid then stops below it. To run the
  grid up to the asymptotic size, give the species parameter data frame a
  `w_max` column equal to `w_inf` as well.
