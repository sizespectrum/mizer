---
name: calibrate-model
description: >-
  Bring a mizer model to steady state and calibrate it to observed data. Use
  whenever the user wants to find the steady state (steady, projectToSteady,
  steadySingleSpecies, steadyNewton), match modelled biomass, numbers, yield or
  growth to observations (calibrateBiomass, matchBiomasses, calibrateNumber,
  matchNumbers, matchGrowth), supply those observations (the
  biomass_observed/biomass_cutoff or number_observed species-parameter columns,
  the yield_observed gear-parameter column), set the level of density-dependent
  reproduction (reproduction_level<-), check convergence (isSteady, getSteadyResidual), or
  diagnose a model that collapses, explodes or will not settle. To ask whether
  the steady state you found is dynamically stable, see the analyse-stability
  skill.
---

# Reaching steady state and calibrating

This covers the tune-to-data loop for an existing `MizerParams` object. To build
the object from scratch first, see the `build-model` skill.

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
| `isSteady(params)` | *(experimental)* ask whether a model **is** at its steady state (boolean) |
| `getSteadyResidual(params)` | *(experimental)* per-capita rates of change across species and sizes, showing where it is not |

During setup and calibration you almost always want `steady()` or
`steadySingleSpecies()`, because holding births constant lets the dynamics settle
reliably onto *a* steady state. Afterwards both re-tune the reproduction
parameters so that density-dependent reproduction reproduces exactly that birth
rate at the new steady state — use the `preserve` argument to choose whether
`reproduction_level` (default), `R_max`, or `erepro` is held fixed during that
re-tuning.

**`steady()` returning is not proof that the model is steady.** Its stopping
test is the relative change in egg production over `t_per`, a proxy for the
drift itself, and the run can also stop on a limit cycle, on a species going
extinct, or simply at `t_max`. All four outcomes return a `MizerParams` object
that looks the same. What distinguishes them is the `"convergence"` attribute,
and the checks in [Verifying the result](#verifying-the-result):

```r
params <- steady(params)
attr(params, "convergence")$type       # "steady", "cycle", "not_converged", "extinction"
attr(params, "convergence")$residual   # largest biomass drift, in 1/year
plot(getSteadyResidual(params))        # which species, and at which sizes
```

`"steady"` with a `residual` above about `0.05` means `tol` was too loose:
tighten it. `"cycle"` means the dynamics settled onto a limit cycle and the
stored state is one point on it — see the `analyse-stability` skill.
`"not_converged"` means the run ran out of time: raise `t_max`, or get a better
starting spectrum from `steadySingleSpecies()` first, and reach for
`steadyNewton()` when the fixed point is dynamically unstable. `"extinction"`
means a species is dying out, which is a modelling problem — see [Diagnosing
calibration problems](#diagnosing-calibration-problems). `steady()` says all of
this in its messages too, but `info_level = 0` suppresses those, so in a script
the attribute is the reliable check.

## The calibration loop

Do this only if you have observations. Observed biomasses live in the
species-parameter column `biomass_observed`, optionally with a
`biomass_cutoff` size threshold below which observations are not counted.
Observed yields live in the gear-parameter column `yield_observed`, which
gives the annual yield of each gear-species pair. The usual order:

```r
params <- steady(params)             # 1. settle onto the steady state
params <- calibrateBiomass(params)   # 2. scale kappa so total modelled biomass
                                     #    matches total observed
params <- matchBiomasses(params)     # 3. adjust each species to its own observation
params <- matchGrowth(params)        # 4. rescale h, gamma, ks, k so each species
                                     #    reaches w_mat / w_inf on schedule
params <- steady(params)             # 5. re-converge after the changes
```

Re-run `steady()` after **any `match…` step** — those rescale each species
separately, which is not a symmetry of the model, so whatever steady state it was
on it is no longer on. They say so when they do it. The `calibrate…()` functions
do **not** need it: they apply one overall scaling factor to the whole model,
which is an exact symmetry and leaves the steady state untouched.

You do not have to remember this. `getSteadyResidual()` answers it, and
`summary(params)` shows the verdict:

```r
summary(params)                        # "biomass drift: 3.2e-05 /year (at steady state)"
plot(getSteadyResidual(params))        # which species, and at which sizes
```

| Function | Adjusts | To match | Breaks steady state |
|---|---|---|---|
| `calibrateBiomass()` | `kappa` (resource level) | total community biomass | no |
| `matchBiomasses()` | per-species abundance | each `biomass_observed` | yes |
| `matchGrowth()` | `h`, `gamma`, `ks`, `k` | von Bertalanffy growth to `w_mat`/`w_inf` | yes |

`matchGrowth()` and `matchBiomasses()` pull on different parameters; alternate
them, re-running `steady()` between, until both are satisfied — usually a few
passes.

**If your observations are numbers rather than weights**, use `calibrateNumber()`
and `matchNumbers()` in place of the two biomass functions. They are the same
functions with the factor of the weight taken out of the size integral, and they
read `number_observed` and `number_cutoff` instead of `biomass_observed` and
`biomass_cutoff`. Mixing the two — matching some species on biomass and others
on numbers — works, because each function ignores the species for which its own
observation column is `NA`.

```r
species_params(params)$number_observed <- c(Cod = 1e6, Herring = 4e8, ...)
params <- calibrateNumber(params)
params <- matchNumbers(params)
params <- steady(params)
```

Whichever pair you use, a `<to>_cutoff` value is the **smallest size the
observation counts**, in grams: a survey that misses fish under 10 g is
`biomass_cutoff = 10`, and the model is then integrated over the same range
rather than over the whole spectrum. Leave it out and the whole size range is
counted, which is the usual reason a model looks like it over-predicts a
species by a wide margin.

**Yields** are not a calibration target in mizer itself: put the observed
annual yield of each gear-species pair into the `yield_observed` column of
`gear_params()` and compare it with the model using
`plotYieldObservedVsModel()`, which sums the observations over the gears. Pass
`gear = ` a subset of the gear names to compare only their catch against only
their observations. Use
`mizerExperimental::matchYield()` to adjust the catchability so that the yields
match. Yields depend on the fishing setup, so make sure gears and effort are
right first — see the `set-up-fishing` skill.

## Density-dependent reproduction

Set how strongly reproduction is density-limited. `reproduction_level` is the
fraction of maximum recruitment realised at steady state (0 = density
independent, closer to 1 = strongly limited):

```r
reproduction_level(params) <- 0.25
```

Alternatively use `setBevertonHolt()` to specify `R_max`, `erepro`, or a per-species
named vector. This does not change the steady state itself — it sets how the
model responds to perturbations away from it. Read the current values back with
`reproduction_level(params)`, useful to check what a model was tuned to before
changing it.

## Verifying the result

```r
isSteady(params)                       # TRUE if settled within tolerance
summary(params)                        # still at the steady state?
plotSpectra(params)                    # sensible, overlapping spectra?
plotGrowthCurves(params, species = "Cod")
plotBiomassObservedVsModel(params)     # points near the 1:1 line?
plotYieldObservedVsModel(params)
```

`project(params, check_steady = TRUE)` makes the same check at the point where it
matters, warning if the run is about to start from a state that is not a fixed
point. It is off by default, because projecting a model away from its steady state
is a normal thing to do.

When the model looks right, project it forward with the `run-simulation` skill
and analyse the results with the `analyse-and-plot` skill.

## Diagnosing calibration problems

- **`steady()` will not settle.** The initial spectrum is probably far from the
  steady state. Run `steadySingleSpecies(params)` first to get a sensible
  starting spectrum, or take a smaller step in whatever parameter you changed and
  re-run. Persistent instability is a case for `steadyNewton()`.
- **A species collapses to near-zero.** Its mortality exceeds the growth it can
  fund. Check its predation-kernel parameters `beta` and `sigma`, its row of the
  interaction matrix, and whether its fishing mortality is too high. See the
  `understand-size-spectrum-dynamics` skill for the underlying physiological
  and trophic mechanics.
- **Biomass matches but growth is wrong (or vice versa).** Alternate
  `matchGrowth()` and `matchBiomasses()`, re-running `steady()` between them.
- **`steady()` said it converged but the results still drift.** Check
  `attr(params, "convergence")$residual`, or `summary(params)`, for how far the
  state actually is from a fixed point, and reduce `steady()`'s `tol` if it is
  too large. `steady()` warns when the two disagree.
- **Results move even though nothing was changed.** The model was not at its
  steady state to begin with. `plot(getSteadyResidual(params))` shows which
  species and which sizes are moving.

## Interactive tuning

For hands-on tuning, the `mizerExperimental` package provides `tuneParams()`, a
Shiny gadget that exposes sliders for the parameters above and re-runs `steady()`
live. It is **not** part of core mizer — install and load `mizerExperimental` to
use it.
