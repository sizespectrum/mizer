# Cheatsheet: Model Setup and Calibration

This cheatsheet gives a quick overview of the functions used to
**build** a mizer model, bring it to **steady state**, **calibrate** it
to observations, and **project** it forward. For full documentation of
each function, follow the links.

Every setter function returns a **new** `MizerParams` object, so always
reassign the result (`params <- steady(params)`). Change species
parameters through
[`species_params(params) <-`](https://sizespectrum.org/mizer/reference/species_params.md)
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
| `gear_params` | fishing gear definitions (see [Cheatsheet: Fishing](https://sizespectrum.org/mizer/articles/cheatsheet-fishing.md)) |
| `no_w`, `min_w`, `max_w` | size-grid resolution and range (`no_w = 100` default) |

Inspect the result with
[`summary(params)`](https://sizespectrum.org/mizer/reference/summary.md),
[`species_params(params)`](https://sizespectrum.org/mizer/reference/species_params.md),
[`getInteraction(params)`](https://sizespectrum.org/mizer/reference/getInteraction.md),
[`gear_params(params)`](https://sizespectrum.org/mizer/reference/gear_params.md)
and
[`gear_params(params)`](https://sizespectrum.org/mizer/reference/resource_params.md).

------------------------------------------------------------------------

## Finding the steady state

A freshly constructed model has only a rough spectrum. Settle it onto a
steady state, which also sets the initial values used by calibration and
[`project()`](https://sizespectrum.org/mizer/reference/project.md).

| Function | Use |
|----|----|
| [`steady(params)`](https://sizespectrum.org/mizer/reference/steady.md) | run the dynamics to convergence with births held fixed (the default) |
| [`projectToSteady(params)`](https://sizespectrum.org/mizer/reference/projectToSteady.md) | version with births responding dynamically; exposes `t_max`, `return_sim` |
| [`steadySingleSpecies(params)`](https://sizespectrum.org/mizer/reference/steadySingleSpecies.md) | set each species to its single-species steady form with births held fixed (fast starting point) |

``` r

params <- steady(params)
```

During model setup and calibration you almost always want
[`steady()`](https://sizespectrum.org/mizer/reference/steady.md) or
[`steadySingleSpecies()`](https://sizespectrum.org/mizer/reference/steadySingleSpecies.md)
because keeping births constant lets the dynamics converge reliably onto
*a* steady state. Afterwards these functions re-tune the reproduction
parameters so that density-dependent reproduction reproduces exactly
that birth rate at the new steady state — use `preserve` to choose
whether `reproduction_level` (default), `R_max`, or `erepro` is held
fixed during that re-tuning.

------------------------------------------------------------------------

## Calibrating to observations

Supply observations in the species-parameter columns `biomass_observed`
and/or `yield_observed`. Then run the calibration loop, **re-running
[`steady()`](https://sizespectrum.org/mizer/reference/steady.md) after
any match/calibrate step**:

| Function | Adjusts | To match |
|----|----|----|
| [`calibrateBiomass(params)`](https://sizespectrum.org/mizer/reference/calibrateBiomass.md) | `kappa` (resource level) | total community biomass |
| [`matchBiomasses(params)`](https://sizespectrum.org/mizer/reference/matchBiomasses.md) | per-species abundance | each `biomass_observed` |
| [`calibrateYield(params)`](https://sizespectrum.org/mizer/reference/calibrateYield.md) | overall abundance scale | total community yield |
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

## When calibration misbehaves

- **A species collapses to near-zero.** Its mortality exceeds the growth
  it can fund. Check its predation-kernel parameters `beta` and `sigma`,
  its row of the interaction matrix, and whether its fishing mortality
  is too high.
- **[`steady()`](https://sizespectrum.org/mizer/reference/steady.md)
  will not settle.** The initial spectrum is probably far from the
  steady state. Run
  [`steadySingleSpecies(params)`](https://sizespectrum.org/mizer/reference/steadySingleSpecies.md)
  first to get a sensible starting spectrum, or take a smaller step in
  whatever parameter you changed and re-run.

Diagnose the fit with

``` r

plotSpectra(params)                    # sensible, overlapping spectra?
plotGrowthCurves(params, species = "Cod")
plotBiomassObservedVsModel(params)     # points near the 1:1 line?
plotYieldObservedVsModel(params)
```

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

## Quick reference

``` r

# ── Build ─────────────────────────────────────────────────────────────────────
params <- newMultispeciesParams(species_params, interaction)
params <- newTraitParams()          # or newCommunityParams(), newSingleSpeciesParams()

# ── Steady state ──────────────────────────────────────────────────────────────
params <- steady(params)
params <- steadySingleSpecies(params)   # fast starting spectrum

# ── Calibrate to data (re-run steady() after each) ────────────────────────────
params <- calibrateBiomass(params)      # total biomass  → kappa
params <- matchBiomasses(params)        # per-species biomass
params <- calibrateYield(params)        # total yield
params <- matchGrowth(params)           # growth → h, gamma, ks, k

# ── Reproduction ──────────────────────────────────────────────────────────────
params <- setBevertonHolt(params, reproduction_level = 0.25)

# ── Project ───────────────────────────────────────────────────────────────────
sim <- project(params, t_max = 20, effort = 1)

# ── Verify ────────────────────────────────────────────────────────────────────
plotBiomassObservedVsModel(params)
plotGrowthCurves(params)
```
