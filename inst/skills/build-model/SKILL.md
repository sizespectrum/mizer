---
name: build-model
description: >-
  Build a new mizer model from a species-parameter data frame. Use
  whenever the user wants to create a MizerParams object with
  newMultispeciesParams() (or newTraitParams, newCommunityParams,
  newSingleSpeciesParams), decide which species-parameter columns to supply and
  which to leave to mizer's allometric defaults, set up an interaction matrix,
  choose the size grid (no_w, min_w, max_w, min_w_pp), or save and reload a
  finished model. Follow this ordered workflow rather than guessing at parameters
  or writing the dynamics by hand. To change a model that already exists see the
  change-parameters skill; fishing is covered by the set-up-fishing skill and
  steady state and calibration by the calibrate-model skill.
---

# Building a mizer model

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
| `w_inf` | von Bertalanffy asymptotic weight (g) — the maximum-size parameter |

Mizer derives defaults for `w_max` (the computational grid boundary, default `1.5 * w_inf`),
`w_repro_max` and `w_mat` from `w_inf`. Everything else has a
sensible default or is calculated.

`w_inf` is the parameter to supply, but it need not be the one you *have*: a
table giving only `w_max`, or only lengths (`l_inf` or `l_max`, with the
length–weight parameters `a` and `b`), is accepted, and mizer fills `w_inf` in
from it and says so. Note that when `w_inf` is taken from `w_max` the two are
equal, so the `1.5 *` headroom above the asymptotic size is not there.

Commonly supplied:

| Column | Meaning |
|---|---|
| `w_min` | Egg size (g, default 0.001) |
| `w_mat` | Maturity weight (g) |
| `beta` | Preferred predator/prey mass ratio (default 30) |
| `sigma` | Width of the lognormal predation kernel (default 2) |
| `k_vb` | von Bertalanffy K — used to derive `h` (and then `gamma`) if `h`/`gamma` absent |
| `h`, `gamma` | Max intake coefficient and search-volume coefficient (alternative to `k_vb`) |
| `a`, `b` | Length–weight conversion parameters ($w = a l^b$, defaults 0.01 and 3) |
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
| `kappa`, `lambda`, `w_pp_cutoff` | resource spectrum coefficient, exponent, and cutoff size |
| `no_w` | number of consumer size bins on the logarithmic grid (default 100) |
| `min_w` | default egg size for species without `w_min` (default 0.001 g); the consumer grid starts at the smallest resulting species egg size |
| `max_w` | largest consumer-grid size; by default the largest species `w_max`, and it cannot be smaller than any species' `w_max` |
| `min_w_pp` | smallest target size for the resource grid; it must be below the consumer grid and resolves to `1e-12` g when omitted |
| `gear_params` | fishing gear definitions (usually omitted and configured later — see the `set-up-fishing` skill; defaults to a knife-edge gear catching every species) |
| `second_order_w` | use the second-order size-advection scheme; see the section "Numerical scheme: watch for numerical diffusion" in the `run-simulation` skill |

The `no_w` bins are equally spaced in log weight between the resulting
consumer-grid minimum and `max_w`. The resource uses the same spacing and adds
smaller bins until `min_w_pp` lies in its smallest bin.

Change gears later with `gear_params(params) <- ...` or `setFishing()` — see the
`set-up-fishing` skill.

## Step 3 — Inspect what you built

```r
summary(params)
species_params(params)     # given + calculated, one row per species
interaction_matrix(params) # the interaction matrix
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

`saveSim()` and `readSim()` do the same for a `MizerSim` object. If the model
needs an extension package, these are the only safe way to persist it — a bare
`readRDS()` silently strips the extension class. See the
`use-extension-packages` skill.

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
