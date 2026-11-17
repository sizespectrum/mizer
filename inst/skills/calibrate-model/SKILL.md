---
name: calibrate-model
description: >-
  Bring a mizer model to steady state and calibrate it to observed data. Use
  whenever the user wants to find the steady state (tuneSteadyState,
  findSteadyState, steadySingleSpecies, and the superseded steady and
  projectToSteady), match modelled biomass, numbers, yield or
  growth to observations (calibrateBiomass, matchBiomasses, calibrateNumber,
  matchNumbers, matchGrowth), supply those observations (the
  biomass_observed/biomass_cutoff or number_observed species-parameter columns,
  the yield_observed gear-parameter column), set the level of density-dependent
  reproduction with the reproduction_level replacement accessor, calculate
  the lifetime reproductive output R0 or ask whether a species can invade
  (getR0, getNPerRecruit), check
  convergence (isSteady, getSteadyResidual), or diagnose a model that collapses,
  explodes or will not settle. To ask whether the steady state you found is
  dynamically stable, see the analyse-stability skill.
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
params <- tuneSteadyState(params)
```

| Function | Use |
|---|---|
| `tuneSteadyState(params)` | solve for the spectra **with births and the resource held fixed**, then adjust the parameters that generate them so those held values are steady too (the default choice while calibrating) |
| `steadySingleSpecies(params)` | set each species to its single-species steady form, births held fixed, without changing the resource — a fast way to get a sensible starting spectrum before `tuneSteadyState()` |
| `findSteadyState(params)` | change **no parameter**: births, resource and spectra settle together wherever the parameters you already have put them |
| `projectUntilSettled(params)` | the same run as `findSteadyState(params)` but returning the `MizerSim`, when you want to watch the approach |
| `isSteady(params)` | *(experimental)* ask whether a model **is** at its steady state (boolean) |
| `getSteadyResidual(params)` | *(experimental)* each size class's contribution to its species' biomass drift, showing where it is not; `rowSums()` gives the per-species drift |

The two finders differ in **what they keep**. `tuneSteadyState()` keeps the
state — the birth rate and the resource abundance you supplied — and moves
`erepro`/`R_max` and the resource capacity `cc_pp` to make it steady.
`findSteadyState()` keeps every parameter and moves the state. During setup and
calibration you almost always want `tuneSteadyState()`, or
`steadySingleSpecies()` first: holding births constant lets the search settle
reliably onto *a* steady state. Use the `preserve` argument to choose whether
`reproduction_level` (default), `R_max`, or `erepro` is held fixed while
`tuneSteadyState()` re-tunes the reproduction parameters.

**When not to use `tuneSteadyState()`.** Holding births fixed is also its
limit: every state it finds has the recruitment the model started with, and it
gets there by changing `erepro`/`R_max` and `cc_pp`. So it is the wrong tool
when

- the question is where a model with *given* parameters settles, for instance
  under a changed fishing effort;
- you are looking for more than one steady state: two alternative stable states
  found with it have different parameters, so they are not states of the same
  model;
- a change in reproduction is the mechanism you are studying, as when predation
  thins a prey species and so raises its recruitment, or when a predator's
  juveniles depend on the small prey that response produces (an emergent Allee
  effect).

Use `findSteadyState()`, or `projectUntilSettled()` to watch the approach, with
the parameters fixed and starting from a suitable initial state, then
`findSteadyState(solver = "newton")` to polish the result. A re-derived `erepro`
greater than 1 is the sign that `tuneSteadyState()` was asked for a state that
does not exist with births held where they were; see
[Diagnosing calibration problems](#diagnosing-calibration-problems).

Both take a `solver` argument. The default `solver = "project"` runs the
dynamics until they settle. `solver = "newton"` *(experimental, needs the
`nleqslv` package)* instead solves the steady-state equation directly, which
converges even when the steady state is dynamically **unstable**, and discovers
the support of the steady state automatically.

**`tuneSteadyState()` returning is not proof that the model is steady.** Its
stopping test is the relative change in egg production over `t_check`, a proxy for
the drift itself, and the run can also stop on a limit cycle, on a species going
extinct, or simply at `t_max`. All four outcomes return a `MizerParams` object
that looks the same. What distinguishes them is the `"convergence"` attribute,
and the checks in [Verifying the result](#verifying-the-result):

```r
params <- tuneSteadyState(params)
attr(params, "convergence")$attractor    # "fixed_point", "limit_cycle" or NA
attr(params, "convergence")$termination  # why the run stopped
attr(params, "convergence")$residual     # largest biomass drift, in 1/year
                                         # (consumers and resource only)
plot(getSteadyResidual(params))          # which species, and at which sizes
```

**`attractor` is the field to read.** It is set from the measured biomass drift
of the model that is returned, so `"fixed_point"` means the drift is within
`residual_tol` (default `0.05`/year) and nothing else does. `termination` says
how the run ended and `converged` whether the solver was happy, neither of which
is a claim about the state.

A `termination` of `"residual_tolerance"` is the ordinary success: the distance
function dropped below `distance_tol` *and* the drift was within `residual_tol`.
`"cycle_detected"` means the search ran into a limit cycle — see the
`analyse-stability` skill. `"time_limit"` means the run ran out of time: raise
`t_max`, or get a better starting spectrum from `steadySingleSpecies()` first,
and reach for `solver = "newton"` when the fixed point is dynamically unstable.
`"extinction"` means a species is dying out, which is a modelling problem — see
[Diagnosing calibration problems](#diagnosing-calibration-problems).
`tuneSteadyState()` says all of this in its messages too, but `info_level = 0`
suppresses those, so in a script the attribute is the reliable check.

A run that hits `t_max` having *met* `distance_tol` is worth recognising: the
message names the drift as the reason, and it means the distance function went
quiet while the model was still moving. Tightening `distance_tol` will not help;
`plot(getSteadyResidual(params))` shows what is actually moving.

## The calibration loop

Do this only if you have observations. Observed biomasses live in the
species-parameter column `biomass_observed`, optionally with a
`biomass_cutoff` size threshold below which observations are not counted.
Observed yields live in the gear-parameter column `yield_observed`, which
gives the annual yield of each gear-species pair. The usual order:

```r
params <- tuneSteadyState(params)    # 1. settle onto the steady state
params <- calibrateBiomass(params)   # 2. scale kappa so total modelled biomass
                                     #    matches total observed
params <- matchBiomasses(params)     # 3. adjust each species to its own observation
params <- matchGrowth(params)        # 4. rescale h, gamma, ks, k so each species
                                     #    reaches w_mat / w_inf on schedule
params <- tuneSteadyState(params)    # 5. re-converge after the changes
```

Re-run `tuneSteadyState()` after **any `match…` step** — those rescale each species
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
them, re-running `tuneSteadyState()` between, until both are satisfied — usually a few
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
params <- tuneSteadyState(params)
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

### Lifetime reproductive output

`getR0(params)` gives the lifetime reproductive output R0 of each species: the
number of eggs that one recruit produces over its life, with growth, mortality,
diffusion and the investment into reproduction held fixed. It ties the
reproduction parameters to the steady state:

- with density-independent reproduction (`R_max = Inf`) every steady state has
  R0 = 1, and because R0 is proportional to `erepro`, dividing `erepro` by R0
  makes the steady state for the current rates a steady state of the model;
- with Beverton-Holt reproduction a steady state has
  R0 = 1 / (1 - `reproduction_level(params)`), not 1.

With Beverton-Holt, Ricker or Shepherd density dependence, or with none, the
recruitment at low abundance equals the egg production, so a rare species
increases if and only if its R0 is above 1. `getR0()` then also answers whether
a species that is absent from a state could invade it. There
`getRDI(params) / getRequiredRDD(params)`, which equals R0 in a steady state,
is 0/0. A custom RDD function whose slope s = dRDD/dRDI at RDI = 0 is not 1
moves the threshold to s R0 > 1, and `getR0()` still returns R0, the egg
output; with `constantRDD()` R0 does not decide whether the species persists. Each rate can be taken from another state, as in
`getR0(params, mort = getMort(other_params))`, to separate the effect of food
from that of predation. `getNPerRecruit(params)` returns the steady-state
abundance per recruit that R0 is calculated from.

Do not calculate R0 by hand from the upwind recursion for the abundance per
recruit. That recursion is right only for the default first-order scheme
without diffusion and is silently off by several per cent otherwise.
`getR0()` solves the same discretised equation as `project()`.

## Verifying the result

```r
isSteady(params)                       # TRUE if settled within tolerance
summary(params)                        # still at the steady state?
attr(getSteadyResidual(params), "other")  # components, if the model has any
plotSpectra(params)                    # sensible, overlapping spectra?
plotGrowthCurves(params, species = "Cod")
plotBiomassObservedVsModel(params)     # points near the 1:1 line?
plotYieldObservedVsModel(params)
```

**A model with components registered by `setComponent()` needs the extra line.**
`isSteady()` and the `summary()` drift cover the consumers and the resource
only: a component's state can be any object at all, so mizer cannot form a
biomass for it, and the functions that find steady states hold it fixed anyway.
Mizer names any component that is moving whenever it reports on the drift, and
`summary()` lists them on their own lines, but the number in front of you is not
about them. `projectUntilSettled()` is what settles them along with everything
else.

When projection holds the size classes above a species' retained support at
zero, the drift diagnostic excludes inflow into those classes if they are
empty. It still includes losses from the last retained class and inflow into
empty classes inside the support. A small biomass drift indicates little net
biomass change in the truncated model; it does not rule out changes within the
size spectrum or establish that `w_max` is large enough for boundary losses to
be negligible.

Populated classes above the support retain their algebraic drift, including
mortality decay. If a projection clears that tail, the instantaneous removal is
not represented by this finite diagnostic. The legacy upwind Euler path without
diffusion can instead fill an empty tail; its diagnostic contributions are also
retained. See the upper-boundary discussion in the numerical details for these
limits.

`project(params, check_steady = TRUE)` makes the same check at the point where it
matters, warning if the run is about to start from a state that is not a fixed
point. It is off by default, because projecting a model away from its steady state
is a normal thing to do.

When the model looks right, project it forward with the `run-simulation` skill
and analyse the results with the `analyse-and-plot` skill.

## Diagnosing calibration problems

- **`tuneSteadyState()` will not settle.** The initial spectrum is probably far from the
  steady state. Run `steadySingleSpecies(params)` first to get a sensible
  starting spectrum, or take a smaller step in whatever parameter you changed and
  re-run. Persistent instability is a case for `solver = "newton"`.
- **A species collapses to near-zero.** Its mortality exceeds the growth it can
  fund. Check its predation-kernel parameters `beta` and `sigma`, its row of the
  interaction matrix, and whether its fishing mortality is too high. See the
  `understand-size-spectrum-dynamics` skill for the underlying physiological
  and trophic mechanics.
- **Biomass matches but growth is wrong (or vice versa).** Alternate
  `matchGrowth()` and `matchBiomasses()`, re-running `tuneSteadyState()`
  between them.
- **`tuneSteadyState()` said it converged but the results still drift.** Check
  `attr(params, "convergence")$residual`, or `summary(params)`, for how far the
  state actually is from a fixed point, and reduce `tuneSteadyState()`'s
  `distance_tol` if it is too large. It warns when the two disagree.
- **`tuneSteadyState()` warns that a species requires "an unrealistic value
  greater than 1 for `erepro`".** The recruitment held during the search needs
  more eggs than that species' spawners can produce in the state that was
  found, and the `"convergence"` attribute will not tell you so: it can still
  report a fixed point. Look at `species_params(params)$erepro`, because the
  warning reads the same for 1.01 as for 1e34. A value far above 1 means the
  state you asked for does not exist with births held at their current values,
  which is what happens when `tuneSteadyState()` is used to look for an
  alternative stable state or for a state that depends on a change in
  reproduction; use `findSteadyState()` or `projectUntilSettled()` with the
  parameters fixed instead (see [Finding the steady
  state](#finding-the-steady-state)). A value just above 1 during calibration
  says the spawners cannot fund the recruitment they have: look at the
  species' feeding level and growth before anything else.
- **Results move even though nothing was changed.** The model was not at its
  steady state to begin with. `plot(getSteadyResidual(params))` shows which
  species and which sizes are moving.

## Interactive tuning

For hands-on tuning, the `mizerExperimental` package provides `tuneParams()`, a
Shiny gadget that exposes sliders for the parameters above and re-runs
`tuneSteadyState()`
live. It is **not** part of core mizer — install and load `mizerExperimental` to
use it.
