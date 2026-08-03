---
name: calibrate-model
description: >-
  Bring a mizer model to steady state and calibrate it to observed data. Use
  whenever the user wants to find the steady state (steady, projectToSteady,
  steadySingleSpecies, steadyNewton), match modelled biomass, yield, or growth to
  observations (calibrateBiomass, matchBiomasses, calibrateYield, matchGrowth),
  set the level of density-dependent reproduction (setBevertonHolt), or diagnose
  why a model will not settle or reproduce observed values.
---

# Bringing a model to steady state and calibrating it

This covers the tune-to-data loop for an existing `MizerParams` object. To build
the object from scratch first, see the `build-multispecies-model` skill.

Every function here **returns a new `MizerParams`** — always reassign
(`params <- f(params, ...)`). Change species parameters through
`species_params(params) <- ...` so dependent quantities recalculate.

## Finding the steady state

A freshly constructed model has only a rough spectrum. Settle it onto a steady
state, which also sets the initial values used by calibration and by `project()`:

```r
params <- steady(params)
```

| Function | Use |
|---|---|
| `steady(params)` | run the dynamics to convergence **with births held fixed**, then store the result as the initial state (the default choice) |
| `steadySingleSpecies(params)` | set each species to its single-species steady form, births held fixed, without changing the resource — a fast way to get a sensible starting spectrum before `steady()` |
| `projectToSteady(params)` | the lower-level routine `steady()` builds on, but with **births responding dynamically**; exposes `t_max`, `tol`, `return_sim` if you need to watch convergence |
| `steadyNewton(params)` | *(experimental)* solve the steady-state equation directly, converging even when the steady state is dynamically unstable |

During setup and calibration you almost always want `steady()` or
`steadySingleSpecies()`, because holding births constant lets the dynamics settle
reliably onto *a* steady state. Afterwards both re-tune the reproduction
parameters so that density-dependent reproduction reproduces exactly that birth
rate at the new steady state — use the `preserve` argument to choose whether
`reproduction_level` (default), `R_max`, or `erepro` is held fixed during that
re-tuning.

## The calibration loop

Do this only if you have observations. They live in the species-parameter
columns `biomass_observed` and/or `yield_observed`, optionally with
`biomass_cutoff` / `yield_cutoff` size thresholds below which observations are
not counted. The usual order:

```r
params <- steady(params)             # 1. settle onto the steady state
params <- calibrateBiomass(params)   # 2. scale kappa so total modelled biomass
                                     #    matches total observed
params <- matchBiomasses(params)     # 3. adjust each species to its own observation
params <- matchGrowth(params)        # 4. rescale h, gamma, ks, k so each species
                                     #    reaches w_mat / w_inf on schedule
params <- steady(params)             # 5. re-converge after the changes
```

Re-run `steady()` after **any** `match…`/`calibrate…` step — those functions move
parameters off the current steady state.

| Function | Adjusts | To match |
|---|---|---|
| `calibrateBiomass()` | `kappa` (resource level) | total community biomass |
| `matchBiomasses()` | per-species abundance | each `biomass_observed` |
| `calibrateYield()` | overall abundance scale | total community yield |
| `matchGrowth()` | `h`, `gamma`, `ks`, `k` | von Bertalanffy growth to `w_mat`/`w_inf` |

`matchGrowth()` and `matchBiomasses()` pull on different parameters; alternate
them, re-running `steady()` between, until both are satisfied — usually a few
passes.

**Yield instead of biomass:** use `calibrateYield()`, which reads
`yield_observed`. Yield calibration depends on the fishing setup, so make sure
gears and effort are right first — see the `set-up-fishing` skill.

## Density-dependent reproduction

Set how strongly reproduction is density-limited. `reproduction_level` is the
fraction of maximum recruitment realised at steady state (0 = density
independent, closer to 1 = strongly limited):

```r
params <- setBevertonHolt(params, reproduction_level = 0.25)
```

Alternatively pass `R_max`, `erepro`, or a per-species named vector. This does
not change the steady state itself — it sets how the model responds to
perturbations away from it. Read the current values back with
`getReproductionLevel(params)`, useful to check what a model was tuned to before
changing it.

## Verifying the result

```r
plotSpectra(params)                    # sensible, overlapping spectra?
plotGrowthCurves(params, species = "Cod")
plotBiomassObservedVsModel(params)     # points near the 1:1 line?
plotYieldObservedVsModel(params)
```

When the model looks right, project it forward with the `run-simulation` skill
and analyse the results with the `analyse-and-plot` skill.

## Diagnosing calibration problems

- **`steady()` will not settle.** The initial spectrum is probably far from the
  steady state. Run `steadySingleSpecies(params)` first to get a sensible
  starting spectrum, or take a smaller step in whatever parameter you changed and
  re-run. Persistent instability is a case for `steadyNewton()`.
- **A species collapses to near-zero.** Its mortality exceeds the growth it can
  fund. Check its predation-kernel parameters `beta` and `sigma`, its row of the
  interaction matrix, and whether its fishing mortality is too high.
- **Biomass matches but growth is wrong (or vice versa).** Alternate
  `matchGrowth()` and `matchBiomasses()`, re-running `steady()` between them.

## Interactive tuning

For hands-on tuning, the `mizerExperimental` package provides `tuneParams()`, a
Shiny gadget that exposes sliders for the parameters above and re-runs `steady()`
live. It is **not** part of core mizer — install and load `mizerExperimental` to
use it.
