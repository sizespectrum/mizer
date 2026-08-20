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

params <- steady(params)
```

| Function | Use |
|----|----|
| [`steady(params)`](https://sizespectrum.org/mizer/reference/steady.md) | run the dynamics to convergence **with births held fixed**, then store the result as the initial state (the default choice) |
| [`steadySingleSpecies(params)`](https://sizespectrum.org/mizer/reference/steadySingleSpecies.md) | set each species to its single-species steady form, births held fixed, without changing the resource — a fast way to get a sensible starting spectrum before [`steady()`](https://sizespectrum.org/mizer/reference/steady.md) |
| [`projectToSteady(params)`](https://sizespectrum.org/mizer/reference/projectToSteady.md) | the lower-level routine [`steady()`](https://sizespectrum.org/mizer/reference/steady.md) builds on, but with **births responding dynamically**; exposes `t_max`, `tol`, `return_sim` if you need to watch convergence |
| [`steadyNewton(params)`](https://sizespectrum.org/mizer/reference/steadyNewton.md) | *(experimental)* solve the steady-state equation directly, converging even when the steady state is dynamically unstable |
| [`isSteady(params)`](https://sizespectrum.org/mizer/reference/isSteady.md) | *(experimental)* ask whether a model **is** at its steady state (boolean) |
| [`getSteadyResidual(params)`](https://sizespectrum.org/mizer/reference/getSteadyResidual.md) | *(experimental)* per-capita rates of change across species and sizes, showing where it is not |

During setup and calibration you almost always want
[`steady()`](https://sizespectrum.org/mizer/reference/steady.md) or
[`steadySingleSpecies()`](https://sizespectrum.org/mizer/reference/steadySingleSpecies.md),
because holding births constant lets the dynamics settle reliably onto
*a* steady state. Afterwards both re-tune the reproduction parameters so
that density-dependent reproduction reproduces exactly that birth rate
at the new steady state — use the `preserve` argument to choose whether
[`reproduction_level`](https://sizespectrum.org/mizer/reference/setBevertonHolt.md)
(default), `R_max`, or `erepro` is held fixed during that re-tuning.

**[`steady()`](https://sizespectrum.org/mizer/reference/steady.md)
returning is not proof that the model is steady.** Its stopping test is
the relative change in egg production over `t_per`, a proxy for the
drift itself, and the run can also stop on a limit cycle, on a species
going extinct, or simply at `t_max`. All four outcomes return a
`MizerParams` object that looks the same. What distinguishes them is the
`"convergence"` attribute, and the checks in [Verifying the
result](#verifying-the-result):

``` r

params <- steady(params)
attr(params, "convergence")$type       # "steady", "cycle", "not_converged", "extinction"
attr(params, "convergence")$residual   # largest biomass drift, in 1/year
plot(getSteadyResidual(params))        # which species, and at which sizes
```

`"steady"` with a `residual` above about `0.05` means `tol` was too
loose: tighten it. `"cycle"` means the dynamics settled onto a limit
cycle and the stored state is one point on it — see the [guide to
analysing dynamic
stability](https://sizespectrum.org/mizer/articles/guide-analyse-stability.md).
`"not_converged"` means the run ran out of time: raise `t_max`, or get a
better starting spectrum from
[`steadySingleSpecies()`](https://sizespectrum.org/mizer/reference/steadySingleSpecies.md)
first, and reach for
[`steadyNewton()`](https://sizespectrum.org/mizer/reference/steadyNewton.md)
when the fixed point is dynamically unstable. `"extinction"` means a
species is dying out, which is a modelling problem — see [Diagnosing
calibration problems](#diagnosing-calibration-problems).
[`steady()`](https://sizespectrum.org/mizer/reference/steady.md) says
all of this in its messages too, but `info_level = 0` suppresses those,
so in a script the attribute is the reliable check.

------------------------------------------------------------------------

## The calibration loop

Do this only if you have observations. Observed biomasses live in the
species-parameter column `biomass_observed`, optionally with a
`biomass_cutoff` size threshold below which observations are not
counted. Observed yields live in the gear-parameter column
`yield_observed`, which gives the annual yield of each gear-species
pair. The usual order:

``` r

params <- steady(params)             # 1. settle onto the steady state
params <- calibrateBiomass(params)   # 2. scale kappa so total modelled biomass
                                     #    matches total observed
params <- matchBiomasses(params)     # 3. adjust each species to its own observation
params <- matchGrowth(params)        # 4. rescale h, gamma, ks, k so each species
                                     #    reaches w_mat / w_inf on schedule
params <- steady(params)             # 5. re-converge after the changes
```

Re-run [`steady()`](https://sizespectrum.org/mizer/reference/steady.md)
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

summary(params)                        # "biomass drift: 3.2e-05 /year (at steady state)"
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
[`steady()`](https://sizespectrum.org/mizer/reference/steady.md)
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
params <- steady(params)
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
plotSpectra(params)                    # sensible, overlapping spectra?
plotGrowthCurves(params, species = "Cod")
plotBiomassObservedVsModel(params)     # points near the 1:1 line?
plotYieldObservedVsModel(params)
```

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

- **[`steady()`](https://sizespectrum.org/mizer/reference/steady.md)
  will not settle.** The initial spectrum is probably far from the
  steady state. Run `steadySingleSpecies(params)` first to get a
  sensible starting spectrum, or take a smaller step in whatever
  parameter you changed and re-run. Persistent instability is a case for
  [`steadyNewton()`](https://sizespectrum.org/mizer/reference/steadyNewton.md).
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
  [`steady()`](https://sizespectrum.org/mizer/reference/steady.md)
  between them.
- **[`steady()`](https://sizespectrum.org/mizer/reference/steady.md)
  said it converged but the results still drift.** Check
  `attr(params, "convergence")$residual`, or `summary(params)`, for how
  far the state actually is from a fixed point, and reduce
  [`steady()`](https://sizespectrum.org/mizer/reference/steady.md)’s
  `tol` if it is too large.
  [`steady()`](https://sizespectrum.org/mizer/reference/steady.md) warns
  when the two disagree.
- **Results move even though nothing was changed.** The model was not at
  its steady state to begin with.
  [`plot(getSteadyResidual(params))`](https://sizespectrum.org/mizer/reference/plot.md)
  shows which species and which sizes are moving.

------------------------------------------------------------------------

## Interactive tuning

For hands-on tuning, the `mizerExperimental` package provides
`tuneParams()`, a Shiny gadget that exposes sliders for the parameters
above and re-runs
[`steady()`](https://sizespectrum.org/mizer/reference/steady.md) live.
It is **not** part of core mizer — install and load `mizerExperimental`
to use it.

------------------------------------------------------------------------

## Quick reference

``` r

# ── Steady state ──────────────────────────────────────────────────────────────
params <- steady(params)
params <- steadySingleSpecies(params)   # fast starting spectrum
params <- steadyNewton(params)          # direct solve (experimental)

# ── Calibrate to data (re-run steady() after each) ────────────────────────────
params <- calibrateBiomass(params)      # total biomass  → kappa
params <- matchBiomasses(params)        # per-species biomass
params <- calibrateNumber(params)       # same, for `number_observed` instead
params <- matchNumbers(params)          #   of `biomass_observed`
params <- matchGrowth(params)           # growth → h, gamma, ks, k
params <- steady(params)                # re-converge

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
