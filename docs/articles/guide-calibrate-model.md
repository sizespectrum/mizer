# Guide: Reaching steady state and calibrating

This guide gives an overview of the functions used to bring a mizer
model to steady state and calibrate it to observed biomasses, yields and
growth. For full documentation of each function, follow the links.

This covers the tune-to-data loop for an existing
[`MizerParams`](https://sizespectrum.org/mizer/reference/MizerParams.md)
object. To build the object from scratch first, see the [guide to
building a mizer
model](https://sizespectrum.org/mizer/articles/guide-build-model.md).

Every function here **returns a new `MizerParams`** — always reassign
(`params <- f(params, ...)`). Change species parameters through
[`species_params(params) <- ...`](https://sizespectrum.org/mizer/reference/species_params.md)
so dependent quantities recalculate.

------------------------------------------------------------------------

## Finding the steady state

A freshly constructed model has only a rough spectrum. Settle it onto a
steady state, which also sets the initial values used by calibration and
by [`project()`](https://sizespectrum.org/mizer/reference/project.md):

``` r

params <- tuneSteadyState(params)
```

| Function | Use |
|----|----|
| [`tuneSteadyState(params)`](https://sizespectrum.org/mizer/reference/tuneSteadyState.md) | solve for the spectra **with births and the resource held fixed**, then adjust the parameters that generate them so those held values are steady too (the default choice while calibrating) |
| [`steadySingleSpecies(params)`](https://sizespectrum.org/mizer/reference/steadySingleSpecies.md) | set each species to its single-species steady form, births held fixed, without changing the resource — a fast way to get a sensible starting spectrum before [`tuneSteadyState()`](https://sizespectrum.org/mizer/reference/tuneSteadyState.md) |
| [`findSteadyState(params)`](https://sizespectrum.org/mizer/reference/findSteadyState.md) | change **no parameter**: births, resource and spectra settle together wherever the parameters you already have put them |
| [`projectUntilSettled(params)`](https://sizespectrum.org/mizer/reference/projectUntilSettled.md) | the same run as [`findSteadyState(params)`](https://sizespectrum.org/mizer/reference/findSteadyState.md) but returning the [`MizerSim`](https://sizespectrum.org/mizer/reference/MizerSim.md), when you want to watch the approach |
| [`isSteady(params)`](https://sizespectrum.org/mizer/reference/isSteady.md) | *(experimental)* `TRUE` if every required measure passes, `FALSE` if any fails, `NA` if coverage is incomplete without a demonstrated failure |
| [`getSteadyAssessment(params)`](https://sizespectrum.org/mizer/reference/getSteadyAssessment.md) | *(experimental)* measures, tolerances, scope and missing support behind that verdict |
| [`getSteadyResidual(params)`](https://sizespectrum.org/mizer/reference/getSteadyResidual.md) | *(experimental)* each size class’s contribution to its species’ biomass drift, showing where it is not; [`rowSums()`](https://rdrr.io/r/base/colSums.html) gives the per-species drift |

The two finders differ in **what they keep**.
[`tuneSteadyState()`](https://sizespectrum.org/mizer/reference/tuneSteadyState.md)
keeps the state — the birth rate and the resource abundance you supplied
— and moves `erepro`/`R_max` and the resource capacity `cc_pp` to make
it steady.
[`findSteadyState()`](https://sizespectrum.org/mizer/reference/findSteadyState.md)
keeps every parameter and moves the state. During setup and calibration
you almost always want
[`tuneSteadyState()`](https://sizespectrum.org/mizer/reference/tuneSteadyState.md),
or
[`steadySingleSpecies()`](https://sizespectrum.org/mizer/reference/steadySingleSpecies.md)
first: holding births constant lets the search settle reliably onto *a*
steady state. Each process supplies a named default tuning strategy.
Built-in Beverton-Holt and Ricker strategies preserve their law-specific
[`reproduction_level`](https://sizespectrum.org/mizer/reference/setBevertonHolt.md).
Use `preserve` for an explicitly supported alternative, such as
`"erepro"` or `"R_max"` for Beverton-Holt, or select named strategies
with `strategies = c(reproduction = "preserve_erepro")`.

Extensions declare which states stay fixed, which evolve, and which
parameters may be adjusted afterwards. A missing declaration is an error
before the search, not permission to guess. `joint` selects registered
strategies that own several coupled processes; those replace the
corresponding local choices. Inspect
[`getProcess(params)`](https://sizespectrum.org/mizer/reference/setProcess.md)
and
[`getComponent(params, component)$definition`](https://sizespectrum.org/mizer/reference/setComponent.md)
for stored declarations. The current Newton solver can hold and later
balance extension states, but cannot evolve arbitrary component states.

Balancing preserves the target state and **all effective rates entering
its equations**, including final assembled recruitment. The independent
checks also enforce declared parameter bounds and preserved quantities.
A failure returns an error, without a generic optimisation fallback or a
second steady-state search. The resulting off-target response, and hence
stability, can still change.

[`balanceRDD(params)`](https://sizespectrum.org/mizer/reference/balanceRDD.md)
preserves the model’s current effective recruitment while adjusting its
current reproduction law. Use `rdd_target` only when a separate
calculation has established the intended constrained target; named
partial targets leave other species’ recruitment and parameters alone
when the strategy supports partial requests.
[`matchGrowth()`](https://sizespectrum.org/mizer/reference/matchGrowth.md),
[`matchBiomasses()`](https://sizespectrum.org/mizer/reference/matchBiomasses.md)
and
[`matchNumbers()`](https://sizespectrum.org/mizer/reference/matchNumbers.md)
retain the chosen law and balance only the selected species.
[`setBevertonHolt()`](https://sizespectrum.org/mizer/reference/setBevertonHolt.md)
is instead an explicit switch to Beverton-Holt, preserving effective
recruitment at the current state. Infeasible explicit controls are
errors rather than silently adjusted values.

**When not to use
[`tuneSteadyState()`](https://sizespectrum.org/mizer/reference/tuneSteadyState.md).**
Holding births fixed is also its limit: every state it finds has the
recruitment the model started with, and it gets there by changing
`erepro`/`R_max` and `cc_pp`. So it is the wrong tool when

- the question is where a model with *given* parameters settles, for
  instance under a changed fishing effort;
- you are looking for more than one steady state: two alternative stable
  states found with it have different parameters, so they are not states
  of the same model;
- a change in reproduction is the mechanism you are studying, as when
  predation thins a prey species and so raises its recruitment, or when
  a predator’s juveniles depend on the small prey that response produces
  (an emergent Allee effect).

Use
[`findSteadyState()`](https://sizespectrum.org/mizer/reference/findSteadyState.md),
or
[`projectUntilSettled()`](https://sizespectrum.org/mizer/reference/projectUntilSettled.md)
to watch the approach, with the parameters fixed and starting from a
suitable initial state, then `findSteadyState(solver = "newton")` to
polish the result. A re-derived `erepro` greater than 1 is the sign that
[`tuneSteadyState()`](https://sizespectrum.org/mizer/reference/tuneSteadyState.md)
was asked for a state that does not exist with births held where they
were; see [Diagnosing calibration
problems](#diagnosing-calibration-problems).

Both take a `solver` argument. The default `solver = "project"` runs the
dynamics until they settle. `solver = "newton"` *(experimental, needs
the `nleqslv` package)* instead solves the steady-state equation
directly, which converges even when the steady state is dynamically
**unstable**, and discovers the support of the steady state
automatically.

**[`tuneSteadyState()`](https://sizespectrum.org/mizer/reference/tuneSteadyState.md)
returning is not proof that the model is steady.** Its stopping test is
the relative change in egg production over `t_check`, a proxy for the
drift itself, and the run can also stop on a limit cycle, on a species
going extinct, or simply at `t_max`. All four outcomes return a
`MizerParams` object that looks the same. What distinguishes them is the
`"convergence"` attribute, and the checks in [Verifying the
result](#verifying-the-result):

``` r

params <- tuneSteadyState(params)
attr(params, "convergence")$attractor    # "fixed_point", "limit_cycle" or NA
attr(params, "convergence")$termination  # why the run stopped
attr(params, "convergence")$assessment   # measures, coverage and scope
attr(params, "convergence")$residual     # separate biomass-drift diagnostic
plot(getSteadyResidual(params))          # which species, and at which sizes
```

**Read `attractor` together with `assessment`.** A `"fixed_point"`
requires all required assessments to pass. With `scope = "conditional"`
that claim applies only with the named states fixed; the report retains
those names. `termination` says how the run ended and `converged`
whether the solver was happy. Neither alone establishes a whole-model
steady state, and the legacy `residual` biomass number is no longer the
verdict.

A `termination` of `"residual_tolerance"` is the ordinary success: the
distance function dropped below `distance_tol` *and* every required
assessment passed. `"cycle_detected"` means the search ran into a limit
cycle — see the [guide to analysing dynamic
stability](https://sizespectrum.org/mizer/articles/guide-analyse-stability.md).
`"time_limit"` means the run ran out of time: raise `t_max`, or get a
better starting spectrum from
[`steadySingleSpecies()`](https://sizespectrum.org/mizer/reference/steadySingleSpecies.md)
first, and reach for `solver = "newton"` when the fixed point is
dynamically unstable. `"extinction"` means a species is dying out, which
is a modelling problem — see [Diagnosing calibration
problems](#diagnosing-calibration-problems).
[`tuneSteadyState()`](https://sizespectrum.org/mizer/reference/tuneSteadyState.md)
says all of this in its messages too, but `info_level = 0` suppresses
those, so in a script the attribute is the reliable check.

A run that hits `t_max` having *met* `distance_tol` is worth
recognising: the message names the drift as the reason, and it means the
distance function went quiet while the model was still moving.
Tightening `distance_tol` will not help;
[`plot(getSteadyResidual(params))`](https://sizespectrum.org/mizer/reference/plot.md)
shows what is actually moving.

------------------------------------------------------------------------

## The calibration loop

Do this only if you have observations. Observed biomasses live in the
species-parameter column `biomass_observed`, optionally with a
`biomass_cutoff` size threshold below which observations are not
counted. Observed yields live in the gear-parameter column
`yield_observed`, which gives the annual yield of each gear-species
pair. The usual order:

``` r

params <- tuneSteadyState(params)    # 1. settle onto the steady state
params <- calibrateBiomass(params)   # 2. scale kappa so total modelled biomass
                                     #    matches total observed
params <- matchBiomasses(params)     # 3. adjust each species to its own observation
params <- matchGrowth(params)        # 4. rescale h, gamma, ks, k so each species
                                     #    reaches w_mat / w_inf on schedule
params <- tuneSteadyState(params)    # 5. re-converge after the changes
```

Re-run
[`tuneSteadyState()`](https://sizespectrum.org/mizer/reference/tuneSteadyState.md)
after **any `match…` step** — those rescale each species separately,
which is not a symmetry of the model, so whatever steady state it was on
it is no longer on. They say so when they do it. The `calibrate…()`
functions do **not** need it: they apply one overall scaling factor to
the whole model, which is an exact symmetry and leaves the steady state
untouched.

You do not have to remember this.
[`getSteadyResidual()`](https://sizespectrum.org/mizer/reference/getSteadyResidual.md)
answers it, and
[`summary(params)`](https://sizespectrum.org/mizer/reference/summary.md)
shows the verdict:

``` r

summary(params)                        # verdict, coverage and tolerance ratios
plot(getSteadyResidual(params))        # which species, and at which sizes
```

| Function | Adjusts | To match | Breaks steady state |
|----|----|----|----|
| [`calibrateBiomass()`](https://sizespectrum.org/mizer/reference/calibrateBiomass.md) | `kappa` (resource level) | total community biomass | no |
| [`matchBiomasses()`](https://sizespectrum.org/mizer/reference/matchBiomasses.md) | per-species abundance | each `biomass_observed` | yes |
| [`matchGrowth()`](https://sizespectrum.org/mizer/reference/matchGrowth.md) | `h`, `gamma`, `ks`, `k` | von Bertalanffy growth to `w_mat`/`w_inf` | yes |

[`matchGrowth()`](https://sizespectrum.org/mizer/reference/matchGrowth.md)
and
[`matchBiomasses()`](https://sizespectrum.org/mizer/reference/matchBiomasses.md)
pull on different parameters; alternate them, re-running
[`tuneSteadyState()`](https://sizespectrum.org/mizer/reference/tuneSteadyState.md)
between, until both are satisfied — usually a few passes.

**If your observations are numbers rather than weights**, use
[`calibrateNumber()`](https://sizespectrum.org/mizer/reference/calibrateNumber.md)
and
[`matchNumbers()`](https://sizespectrum.org/mizer/reference/matchNumbers.md)
in place of the two biomass functions. They are the same functions with
the factor of the weight taken out of the size integral, and they read
`number_observed` and `number_cutoff` instead of `biomass_observed` and
`biomass_cutoff`. Mixing the two — matching some species on biomass and
others on numbers — works, because each function ignores the species for
which its own observation column is `NA`.

``` r

species_params(params)$number_observed <- c(Cod = 1e6, Herring = 4e8, ...)
params <- calibrateNumber(params)
params <- matchNumbers(params)
params <- tuneSteadyState(params)
```

Whichever pair you use, a `<to>_cutoff` value is the **smallest size the
observation counts**, in grams: a survey that misses fish under 10 g is
`biomass_cutoff = 10`, and the model is then integrated over the same
range rather than over the whole spectrum. Leave it out and the whole
size range is counted, which is the usual reason a model looks like it
over-predicts a species by a wide margin.

**Yields** are not a calibration target in mizer itself: put the
observed annual yield of each gear-species pair into the
`yield_observed` column of
[`gear_params()`](https://sizespectrum.org/mizer/reference/gear_params.md)
and compare it with the model using
[`plotYieldObservedVsModel()`](https://sizespectrum.org/mizer/reference/plotYieldObservedVsModel.md),
which sums the observations over the gears. Pass `gear =` a subset of
the gear names to compare only their catch against only their
observations. Use
[`mizerExperimental::matchYield()`](https://sizespectrum.org/mizerExperimental/reference/matchYield.html)
to adjust the catchability so that the yields match. Yields depend on
the fishing setup, so make sure gears and effort are right first — see
the [guide to setting up
fishing](https://sizespectrum.org/mizer/articles/guide-set-up-fishing.md).

------------------------------------------------------------------------

## Density-dependent reproduction

Set how strongly reproduction is density-limited. `reproduction_level`
is the fraction of maximum recruitment realised at steady state (0 =
density independent, closer to 1 = strongly limited):

``` r

reproduction_level(params) <- 0.25
```

Alternatively use
[`setBevertonHolt()`](https://sizespectrum.org/mizer/reference/setBevertonHolt.md)
to specify `R_max`, `erepro`, or a per-species named vector. This does
not change the steady state itself — it sets how the model responds to
perturbations away from it. Read the current values back with
`reproduction_level(params)`, useful to check what a model was tuned to
before changing it.

------------------------------------------------------------------------

## Verifying the result

``` r

isSteady(params)                       # TRUE if settled within tolerance
summary(params)                        # still at the steady state?
attr(getSteadyResidual(params), "other")  # components, if the model has any
plotSpectra(params)                    # sensible, overlapping spectra?
plotGrowthCurves(params, species = "Cod")
plotBiomassObservedVsModel(params)     # points near the 1:1 line?
plotYieldObservedVsModel(params)
```

Whole-model coverage is the default. For built-in abundance states, the
default measure integrates the magnitude of each size-class derivative,
weighted by body mass. Each species and the resource is checked
separately against `atol + tol * integral(w * abs(N))`, with
`atol = 1e-12` g/year and `tol = 0.05` per year by default. The model’s
size quadrature is respected. At zero abundance, zero change passes and
positive inflow must fit the absolute tolerance. Taking magnitudes
before summing prevents changes at different sizes from cancelling.

For the previous net-biomass criterion, explicitly request
`getSteadyAssessment(params, criterion = "biomass")`. This can pass
while the shape of a spectrum is changing, so it is a deliberate weaker
criterion. Extension measures have their own units and tolerances; `tol`
applies to the built-in measures. Use named `tolerances` overrides keyed
by `owner/measure` to adjust an individual reported measure.

A component needs its own convergence capability even if it supplies an
equation or a biomass measurement. Missing support appears in
`assessment$missing` and `assessment$unavailable`; it cannot be treated
as passing. Callback errors are errors. Solvers stop before starting if
required support is missing.

Use `scope = "conditional", fixed = "component:pool"` to analyse a
subsystem with a particular dynamic component held fixed. Other IDs are
`consumers` and `resource`; several can be supplied together. The
current solvers cannot hold consumers fixed, and Newton cannot evolve
arbitrary extension states. A conditional result remains labelled
conditional. A structurally prescribed input is recorded separately from
a dynamic state fixed for this analysis.

When projection holds the size classes above a species’ retained support
at zero, the drift diagnostic excludes inflow into those classes if they
are empty. It still includes losses from the last retained class and
inflow into empty classes inside the support. A small biomass drift
indicates little net biomass change in the truncated model; it does not
rule out changes within the size spectrum or establish that `w_max` is
large enough for boundary losses to be negligible.

Populated classes above the support retain their algebraic drift,
including mortality decay. If a projection clears that tail, the
instantaneous removal is not represented by this finite diagnostic. The
legacy upwind Euler path without diffusion can instead fill an empty
tail; its diagnostic contributions are also retained. See the
upper-boundary discussion in the numerical details for these limits.

`project(params, check_steady = TRUE)` makes the same check at the point
where it matters, warning if the run is about to start from a state that
is not a fixed point. It is off by default, because projecting a model
away from its steady state is a normal thing to do.

When the model looks right, project it forward with the [guide to
running a mizer
simulation](https://sizespectrum.org/mizer/articles/guide-run-simulation.md)
and analyse the results with the [guide to analysing and plotting mizer
results](https://sizespectrum.org/mizer/articles/guide-analyse-and-plot.md).

------------------------------------------------------------------------

## Diagnosing calibration problems

- **[`tuneSteadyState()`](https://sizespectrum.org/mizer/reference/tuneSteadyState.md)
  will not settle.** The initial spectrum is probably far from the
  steady state. Run `steadySingleSpecies(params)` first to get a
  sensible starting spectrum, or take a smaller step in whatever
  parameter you changed and re-run. Persistent instability is a case for
  `solver = "newton"`.
- **A species collapses to near-zero.** Its mortality exceeds the growth
  it can fund. Check its predation-kernel parameters `beta` and `sigma`,
  its row of the interaction matrix, and whether its fishing mortality
  is too high. See the [guide to understanding size-spectrum
  dynamics](https://sizespectrum.org/mizer/articles/guide-understand-size-spectrum-dynamics.md)
  for the underlying physiological and trophic mechanics.
- **Biomass matches but growth is wrong (or vice versa).** Alternate
  [`matchGrowth()`](https://sizespectrum.org/mizer/reference/matchGrowth.md)
  and
  [`matchBiomasses()`](https://sizespectrum.org/mizer/reference/matchBiomasses.md),
  re-running
  [`tuneSteadyState()`](https://sizespectrum.org/mizer/reference/tuneSteadyState.md)
  between them.
- **[`tuneSteadyState()`](https://sizespectrum.org/mizer/reference/tuneSteadyState.md)
  said it converged but the results still drift.** Check
  `attr(params, "convergence")$assessment`, or `summary(params)`, for
  which measures still fail, and reduce
  [`tuneSteadyState()`](https://sizespectrum.org/mizer/reference/tuneSteadyState.md)’s
  `distance_tol` if it is too large. It warns when the two disagree.
- **[`tuneSteadyState()`](https://sizespectrum.org/mizer/reference/tuneSteadyState.md)
  warns that a species requires “an unrealistic value greater than 1 for
  `erepro`”.** The recruitment held during the search needs more eggs
  than that species’ spawners can produce in the state that was found,
  and the `"convergence"` attribute will not tell you so: it can still
  report a fixed point. Look at `species_params(params)$erepro`, because
  the warning reads the same for 1.01 as for 1e34. A value far above 1
  means the state you asked for does not exist with births held at their
  current values, which is what happens when
  [`tuneSteadyState()`](https://sizespectrum.org/mizer/reference/tuneSteadyState.md)
  is used to look for an alternative stable state or for a state that
  depends on a change in reproduction; use
  [`findSteadyState()`](https://sizespectrum.org/mizer/reference/findSteadyState.md)
  or
  [`projectUntilSettled()`](https://sizespectrum.org/mizer/reference/projectUntilSettled.md)
  with the parameters fixed instead (see [Finding the steady
  state](#finding-the-steady-state)). A value just above 1 during
  calibration says the spawners cannot fund the recruitment they have:
  look at the species’ feeding level and growth before anything else.
- **Results move even though nothing was changed.** The model was not at
  its steady state to begin with. `plot(getSteadyResidual(params))`
  shows which species and which sizes are moving.

------------------------------------------------------------------------

## Interactive tuning

For hands-on tuning, the `mizerExperimental` package provides
`tuneParams()`, a Shiny gadget that exposes sliders for the parameters
above and re-runs
[`tuneSteadyState()`](https://sizespectrum.org/mizer/reference/tuneSteadyState.md)
live. It is **not** part of core mizer — install and load
`mizerExperimental` to use it.

------------------------------------------------------------------------

## Quick reference

``` r

# ── Steady state ──────────────────────────────────────────────────────────────
params <- tuneSteadyState(params)
params <- steadySingleSpecies(params)   # fast starting spectrum
params <- tuneSteadyState(params, solver = "newton")  # direct solve (experimental)

# ── Calibrate to data (re-run tuneSteadyState() after each) ───────────────────
params <- calibrateBiomass(params)      # total biomass  → kappa
params <- matchBiomasses(params)        # per-species biomass
params <- calibrateNumber(params)       # same, for `number_observed` instead
params <- matchNumbers(params)          #   of `biomass_observed`
params <- matchGrowth(params)           # growth → h, gamma, ks, k
params <- tuneSteadyState(params)       # re-converge

# ── Reproduction ──────────────────────────────────────────────────────────────
reproduction_level(params) <- 0.25
reproduction_level(params)              # what the model is currently tuned to

# ── Verify ────────────────────────────────────────────────────────────────────
isSteady(params)                        # TRUE if settled within tolerance
summary(params)                         # includes the biomass-drift verdict
plot(getSteadyResidual(params))         # which species and sizes are still moving
plotSpectra(params)
plotGrowthCurves(params)
plotBiomassObservedVsModel(params)
plotYieldObservedVsModel(params)
```
