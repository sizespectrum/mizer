# Cheatsheet: Model Setup and Calibration

This cheatsheet gives a quick overview of the functions used to
**build** a mizer model, bring it to **steady state**, **calibrate** it
to observations, and **project** it forward. For full documentation of
each function, follow the links.

Every setter function returns a **new** `MizerParams` object, so always
reassign the result (`params <- steady(params)`). Change species
parameters through
[`given_species_params(params) <-`](https://sizespectrum.org/mizer/reference/species_params.md)
so that dependent quantities are recalculated.

------------------------------------------------------------------------

## Creating a model

Choose a constructor for the type of model you need.

| Function | Model type |
|----|----|
| [`newSingleSpeciesParams()`](https://sizespectrum.org/mizer/reference/newSingleSpeciesParams.md) | one species in a fixed background |
| [`newCommunityParams()`](https://sizespectrum.org/mizer/reference/newCommunityParams.md) | single size spectrum, no species identity |
| [`newTraitParams()`](https://sizespectrum.org/mizer/reference/newTraitParams.md) | several species differing only in asymptotic size |
| [`newMultispeciesParams()`](https://sizespectrum.org/mizer/reference/newMultispeciesParams.md) | fully general multi-species model |

Most work uses
[`newMultispeciesParams()`](https://sizespectrum.org/mizer/reference/newMultispeciesParams.md),
driven by a **species parameter data frame** — one row per species. Only
two columns are required:

| Column | Meaning |
|----|----|
| `species` | species name |
| `w_inf` | von Bertalanffy asymptotic weight (g) — the required maximum-size parameter |

`w_max` (the computational grid boundary) defaults to `1.5 * w_inf`;
`w_mat`, `beta`, `sigma`, `h`, `gamma`, `alpha`, `erepro`, `R_max`, …
all have defaults or are calculated. Weights are in **grams**, lengths
in **cm**, time in **years**.

``` r

species_params <- read.csv(
    system.file("extdata", "NS_species_params.csv", package = "mizer"))
params <- newMultispeciesParams(species_params)
```

Useful optional arguments to
[`newMultispeciesParams()`](https://sizespectrum.org/mizer/reference/newMultispeciesParams.md):

| Argument | Effect |
|----|----|
| `interaction` | species × species matrix of overlaps in `[0, 1]` (default all 1) |
| `gear_params` | fishing gear definitions (see the fishing cheatsheet) |
| `no_w`, `min_w`, `max_w` | size-grid resolution and range (`no_w = 100` default) |

Inspect the result with
[`summary(params)`](https://sizespectrum.org/mizer/reference/summary.md),
[`species_params(params)`](https://sizespectrum.org/mizer/reference/species_params.md),
[`getInteraction(params)`](https://sizespectrum.org/mizer/reference/getInteraction.md),
and
[`gear_params(params)`](https://sizespectrum.org/mizer/reference/gear_params.md).

------------------------------------------------------------------------

## Finding the steady state

A freshly constructed model has only a rough spectrum. Settle it onto a
steady state, which also sets the initial values used by calibration and
[`project()`](https://sizespectrum.org/mizer/reference/project.md).

| Function | Use |
|----|----|
| [`steady(params)`](https://sizespectrum.org/mizer/reference/steady.md) | run the dynamics to convergence (the default) |
| [`projectToSteady(params)`](https://sizespectrum.org/mizer/reference/projectToSteady.md) | lower-level version exposing `t_max`, `tol`, `return_sim` |
| [`steadySingleSpecies(params)`](https://sizespectrum.org/mizer/reference/steadySingleSpecies.md) | set each species to its single-species steady form (fast starting point) |
| [`steadyNewton(params)`](https://sizespectrum.org/mizer/reference/steadyNewton.md) | solve the steady-state equation directly; converges even when unstable |

``` r

params <- steady(params)
```

------------------------------------------------------------------------

## Calibrating to observations

Supply observations in the species-parameter columns `biomass_observed`
and/or `yield_observed` (optionally with `biomass_cutoff` /
`yield_cutoff` size thresholds). Then run the calibration loop,
**re-running
[`steady()`](https://sizespectrum.org/mizer/reference/steady.md) after
any match/calibrate step**:

| Function | Adjusts | To match |
|----|----|----|
| [`calibrateBiomass(params)`](https://sizespectrum.org/mizer/reference/calibrateBiomass.md) | `kappa` (resource level) | total community biomass |
| [`matchBiomasses(params)`](https://sizespectrum.org/mizer/reference/matchBiomasses.md) | per-species abundance | each `biomass_observed` |
| [`calibrateYield(params)`](https://sizespectrum.org/mizer/reference/calibrateYield.md) | overall abundance scale | total community yield |
| [`matchYields(params)`](https://sizespectrum.org/mizer/reference/matchYields.md) | per-species abundance | each `yield_observed` |
| [`matchGrowth(params)`](https://sizespectrum.org/mizer/reference/matchGrowth.md) | `h`, `gamma`, `ks`, `k` | growth to `w_mat` / `w_inf` |

``` r

params <- calibrateBiomass(params)   # total biomass
params <- matchBiomasses(params)     # per species
params <- matchGrowth(params)        # growth
params <- steady(params)             # re-converge
```

[`matchGrowth()`](https://sizespectrum.org/mizer/reference/matchGrowth.md)
and
[`matchBiomasses()`](https://sizespectrum.org/mizer/reference/matchBiomasses.md)
pull on different parameters; alternate them, re-running
[`steady()`](https://sizespectrum.org/mizer/reference/steady.md)
between, until both are satisfied.

------------------------------------------------------------------------

## Density-dependent reproduction

[`setBevertonHolt()`](https://sizespectrum.org/mizer/reference/setBevertonHolt.md)
sets how strongly reproduction is density-limited. `reproduction_level`
is the fraction of maximum recruitment realised at steady state (0 =
density independent, closer to 1 = strongly limited).

``` r

params <- setBevertonHolt(params, reproduction_level = 0.25)
```

You can instead pass `R_max`, `erepro`, or a per-species named vector.

------------------------------------------------------------------------

## Projecting forward

[`project()`](https://sizespectrum.org/mizer/reference/project.md) runs
the model and returns a `MizerSim`. See the analysis-and-plotting
cheatsheet for exploring the result.

``` r

sim <- project(params, t_max = 20, effort = 1)
```

| Argument | Meaning |
|----|----|
| `effort` | scalar, per-gear vector, or time × gear array (see fishing cheatsheet) |
| `t_max` | number of years to simulate |
| `dt` | integration time step (reduce if a run is unstable) |
| `t_save` | interval at which output is stored |
| `method` | `"euler"` (default), `"predictor_corrector"`, `"tr_bdf2"` |

Pass a `MizerSim` back to
[`project()`](https://sizespectrum.org/mizer/reference/project.md) to
continue from where it ended.

------------------------------------------------------------------------

## Verifying the model

``` r

plotSpectra(params)                    # sensible, overlapping spectra?
plotGrowthCurves(params, species = "Cod")
plotBiomassObservedVsModel(params)     # points near the 1:1 line?
plotYieldObservedVsModel(params)
```

------------------------------------------------------------------------

## Quick reference

``` r

# ── Build ─────────────────────────────────────────────────────────────────────
params <- newMultispeciesParams(species_params, interaction)
params <- newTraitParams()          # or newCommunityParams(), newSingleSpeciesParams()

# ── Steady state ──────────────────────────────────────────────────────────────
params <- steady(params)
params <- steadySingleSpecies(params)   # fast starting spectrum
params <- steadyNewton(params)          # direct solve (experimental)

# ── Calibrate to data (re-run steady() after each) ────────────────────────────
params <- calibrateBiomass(params)      # total biomass  → kappa
params <- matchBiomasses(params)        # per-species biomass
params <- calibrateYield(params)        # total yield
params <- matchYields(params)           # per-species yield
params <- matchGrowth(params)           # growth → h, gamma, ks, k

# ── Reproduction ──────────────────────────────────────────────────────────────
params <- setBevertonHolt(params, reproduction_level = 0.25)

# ── Project ───────────────────────────────────────────────────────────────────
sim <- project(params, t_max = 20, effort = 1)

# ── Verify ────────────────────────────────────────────────────────────────────
plotBiomassObservedVsModel(params)
plotGrowthCurves(params)
```
