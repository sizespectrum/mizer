---
name: analyse-stability
description: >-
  Analyse the dynamic stability of a mizer steady state and characterise limit
  cycles (experimental). Use whenever the user asks whether a steady state is
  stable or unstable, wants the leading eigenvalue, the period of an emergent
  oscillation, a Hopf bifurcation, the spectral radius of the numerical time
  step, a limit cycle to build or plot, or a bifurcation diagram over fishing
  effort — via getStability(), getDiscreteStability(), findSteadyState(),
  getOscillationModeSim() and scanModel(). This
  skill and calibrate-model share findSteadyState() and getSteadyResidual(): use
  calibrate-model to find a steady state, this skill to ask whether the state you
  found is stable. Assumes the standard semichemostat resource dynamics.
---

# Analysing dynamic stability

Tools for asking whether a mizer steady state is dynamically **stable** — and,
when it is not, characterising the **limit cycle** that replaces it. These are
**experimental**: their interface may still change. `findSteadyState(solver =
"newton")` assumes the
standard semichemostat resource dynamics and solves for the resource alongside
the fish, so the resource density and the feeding levels it implies are
self-consistent even where consumers are satiated. `getStability()` and
`getDiscreteStability()` include every supported dynamic resource coordinate
alongside the fish. Continuous analysis needs a resource equation capability;
discrete analysis uses the supported resource step. Both require an appropriate
autonomous-model contract for the providers used by the analysis.

Distinct from calibration: the `calibrate-model` skill gets you *onto* a fixed
point; what follows analyses the *dynamics around* it. The natural entry point
is `findSteadyState(params, solver = "newton")`, which converges even when the
fixed point is dynamically unstable, unlike the default `solver = "project"`.

## Is the steady state stable? — `getStability()`

`getStability(params)` computes the continuous-time eigenvalues of the model.
Mizer discretises the size axis but not time, so on the size grid the model is a
system of ODEs; `getStability()` differentiates their right-hand side directly
and takes the eigenvalues of that Jacobian. **No time step is involved**, so the
eigenvalues are a property of the model rather than of a solver. The result
reports whether the state is stable, unstable or inconclusive, the maximum real part of the
eigenvalues (`max_real_part`; > 0 means unstable), and the **period** at which
the model rings (`dominant_period` for the dominant mode,
`oscillation_period` for the leading *oscillatory* one). A complex pair is an
oscillatory mode, not a bifurcation: its period is the period of an emerging
limit cycle only where its real part is zero, which is a Hopf bifurcation and
which a single spectrum cannot show. Watch
`Re(leading_oscillatory_eigenvalue)` cross zero over a `scanModel()` before
calling one. It also returns `leading_eigenvectors`, the
top two eigenvectors of the state space split into `$fish`, a complex array
`(n_species, n_sizes, 2)`, and `$resource`, a complex matrix `(n_w_full, 2)`,
each normalised to maximum modulus 1 (the spatial shape of the oscillation).

```r
params <- findSteadyState(params, solver = "newton")  # sit exactly on the (possibly unstable) fixed point
stab   <- getStability(params)
stab                                     # stable/unstable, growth rate, cycle period
```

- `effort` sets the fishing effort used when forming the Jacobian.
- `h` is the relative finite-difference step. Re-running with a different `h` is
  the cheapest check that the model is smooth enough for the analysis to mean
  anything: if the answer moves, do not trust it.

**Both `getStability()` and `getOscillationModeSim()` use the state stored in
the object.** At a state that fails its scoped steady-state assessment, raw
Jacobian diagnostics do not establish equilibrium stability: the verdict is
inconclusive. `getOscillationModeSim()` requires complete, current support and
an autonomous equilibrium for its recorded scope. Run `findSteadyState()` or
`tuneSteadyState()` first, then inspect `getSteadyAssessment()` — see the
`calibrate-model` skill. An old stability list without coverage and scope must
be recomputed.

`tuneSteadyState()` and `findSteadyState()` attach a related `"convergence"`
attribute, which answers three questions separately: `termination` says why the
run stopped, `converged` whether the solver met its own criterion, and
`attractor` what the state reached actually is — `"fixed_point"`,
`"limit_cycle"` or `NA`. **Read `attractor` together with `assessment`**:
a fixed-point label requires the measures for its recorded scope to pass.
`residual` remains a separate biomass-drift diagnostic. An explicitly conditional
fixed point is not a whole-model equilibrium claim.

Stopping on a fixed point takes both a distance function below `distance_tol`
and every required assessment passing; missing assessment support stops the
operation before the search. Limit cycles are detected from the per-species biomass
series sampled at every time step, at every check rather than only when the
distance criterion has failed, because a cycle whose period divides `t_check` is
sampled at one phase and looks converged. The relative-amplitude floor for calling an
oscillation a cycle is the `amplitude_tol` argument (default `0.01`).

### What the analysis covers

`getSteadyAssessment()` and `isSteady()` cover the whole model by default,
including registered component assessments. A conditional request names the
states held fixed, and its scope stays attached to the result. The current
Newton solver can solve consumers and resource only; dynamic component states
must be explicitly held fixed when using it. The projection solver can advance
components, provided their convergence measures are declared.

Stability means the whole model by default. Dynamic component Jacobians are
not available in this release, so a model containing an uncovered dynamic
component returns an incomplete, inconclusive report. It does not silently
compute a subsystem and label it as the whole model. Request the subsystem
explicitly when that is the question you want to ask:

```r
stab <- getStability(params, scope = "conditional", fixed = "component:pool")
stab$scope
stab$coverage
```

The fixed owners remain recorded throughout the result and any generated
oscillation-mode simulation. A fixed dynamic pool is different from an input
declared by the model. Explicit prescribed-input declarations, including
coordinate masks for partially prescribed states, remove those coordinates
from both continuous and discrete eigenproblems. Built-in `resource_constant`
declares the resource as prescribed. Excluded coordinates are reconstructed as
zero perturbations in eigenvectors and stay unchanged in mode simulations.

Mizer never decides that a coordinate is prescribed because its abundance,
derivative or Jacobian row happens to be zero. Entirely absent fish species still
have invasion modes, and genuine zero eigenvalues or unit multipliers remain in
the spectrum. This differs from the active set used by the Newton solver.

Read `status`, `stable`, `coverage` and `assessment` together:

| Evidence for the requested scope | `status` | `stable` |
|---|---|---|
| Complete autonomous equilibrium, all modes strictly decay | `stable` | `TRUE` |
| Complete autonomous equilibrium, a mode grows | `unstable` | `FALSE` |
| Marginal modes, missing support, or a non-equilibrium state | `inconclusive` | `NA` |
| Complete coverage with no dynamic coordinates | `not_applicable` | `NA` |

`stability_tol` sets the numerical margin around zero real part or unit modulus.
The result records this tolerance, the evaluation time, prescribed coordinates,
conditional exclusions and reasons for missing support. Time-dependent forcing
needs its own suitable analysis; a stationary state alone does not establish
that its Jacobian is autonomous. `isTRUE(stab$stable)` safely tests for an
established stable verdict.

## Is mizer's time step stable? — `getDiscreteStability()`

A different question, and the one to ask when a simulation disagrees with
`getStability()`. `getDiscreteStability(params, dt = 1)` linearises the one-step
map that `project(method = "euler")` takes at that step size and returns its
`discrete_eigenvalues` and `spectral_radius`. Below 1 means the *numerical
scheme* does not amplify perturbations at that `dt`. It is that step and not an
approximation to it: the same transport coefficients, the same Thomas solver and
the model's own `resource_dynamics`.

With excluded coordinates, this is the full step with those inputs reset to
their reference values between steps, followed by selection of the dynamic
outputs. It is not a newly derived numerical solver for the constrained
continuous equations; an implicit step can couple coordinates internally.

The two verdicts can differ, and that is the point: the implicit transport solve
damps oscillations artificially, so an unstable steady state can have a spectral
radius below 1 at a large `dt` and a simulation will then settle onto a state the
model does not hold. If `getStability()` says unstable but the simulation goes
flat, check `getDiscreteStability()` at the `dt` you are projecting with, and
reduce `dt` (or use `method = "second_order"`).

The discrete eigenvalues cannot be converted into continuous ones by any exact
algebraic relation, because mizer's step is not fully implicit: the rates are
evaluated at the state at the start of the step and only the transport solve is
implicit. That is why `getStability()` differentiates the rates of change
directly instead of inverting the one-step map.

## Visualising the oscillation — `getOscillationModeSim()`

`getOscillationModeSim(params)` takes the output of `findSteadyState()` and builds a
`MizerSim` covering **exactly one period** of the leading oscillatory mode in
the linear approximation,

\[ x(t) = x^* + A\,\mathrm{Re}\!\left[e^{i\omega t}\,\mathbf v\right], \]

over the whole state \(x = (N, n_{pp})\), where \(\mathbf v\) is
`getStability()`'s `leading_oscillatory_eigenvector` — the eigenvector of the dominant
*oscillatory* mode, which is not in general the dominant mode — and \(A\) is
scaled so that the largest relative swing in **species biomass** equals the
`amplitude` argument (default 10%): the hardest-oscillating species departs
that far from its steady biomass and no species departs further. Fish and
resource share one \(A\), so the resource oscillates with the phase the mode
gives it; on `NS_params` at effort 1.5 it leads the fish biomass by about 0.3
years of the 5.12-year period. The growth factor \(e^{\sigma t}\) is deliberately
omitted so the oscillation closes after one period; the run ends exactly at
\(T\), with a shortened final interval when `t_save` does not divide it, so the
last state is the first state again. \(\sigma\) is recorded as `growth_rate` in
the result's `sim_params`, and it is what decides whether this is a cycle the
model settles onto or a transient it rings with on the way back.

A biomass amplitude does not bound the individual size classes — a cohort
trough is a far larger relative excursion than the biomass integral over it —
so `getOscillationModeSim()` reports when the perturbation drives an abundance
negative and it has to be clipped at zero. On `NS_params` at effort 1.5 that
starts at around `amplitude = 0.2`. Past that point the picture is no longer
the linear mode. The result is an ordinary `MizerSim`, so plot it with the standard
tools (see the `analyse-and-plot` skill):

```r
params <- findSteadyState(params, solver = "newton")
sim    <- getOscillationModeSim(params, amplitude = 0.1)
plotBiomass(sim)                         # biomass oscillation over one period
animate(plotSpectra(sim))               # the travelling wave in the spectrum
```

This is the *linear* mode — the shape near onset. For the fully nonlinear cycle,
project the real dynamics with the second-order scheme (see the
`run-simulation` skill).

## Where does the cycle appear? — `scanModel()`

`scanModel()` sweeps any aspect of the model and, at each value, follows the
attractor of the full dynamics and records the long-term **range** of a summary
quantity. Draw that range with `style = "envelope"` and you have a bifurcation
diagram: a stable steady state shows as a single line, a limit cycle as a **band**
between its minimum and maximum, so a Hopf bifurcation appears as the value at
which the band opens up.

```r
scan <- scanModel(params, scan_values = seq(0, 1.5, 0.05),
                  set_func = scanEffort(), value_func = getYield)
plot(scan, style = "envelope")
```

- `set_func` says what to vary — `scanEffort()` for the fishing effort of every
  gear together, `scanFishingMortality()` for the fishing mortality on one
  species, `scanSpeciesParam()` for any species parameter. `value_func` says what
  to measure: any of `getBiomass()`, `getYield()`, `getSSB()`, `getN()`. `species`
  restricts which series are shown. See the `analyse-and-plot` skill.
- The settling stage runs `projectUntilSettled()`, whose `distance_tol`,
  `residual_tol`, `amplitude_tol` and `extinction_threshold` are exposed for
  tuning when a run is slow to settle or a species is collapsing.
- The returned `MizerScan` object is a data frame, so the numbers behind the
  diagram are available without re-running the scan.

### Reading and subsetting a scan

The first two columns of a `MizerScan` have **dynamic names**: column 1 is named
by `scan_name` (for example, `"Fishing effort"`), and column 2 by `value_name`
(for example, `"Yield"`). Read them by position in reusable code, or obtain their
names with `names(scan)[1:2]`:

```r
scanned <- scan[[1]]
measured <- scan[[2]]
fixed <- scan[scan$attractor %in% "fixed_point", ]
plot(fixed)
```

The other columns have fixed names: `Species` identifies the series; `ymin`
and `ymax` give its range; `attractor` is `"fixed_point"`, `"limit_cycle"` or
`NA`; `termination` and `converged` describe how the run stopped; and `period`
and `residual` report the cycle period and remaining biomass drift. Use
`attractor`, rather than `converged` alone, to identify equilibrium or cycle
results.

Use `[` for row subsetting. It preserves the plotting metadata and recomputes
`at_max` and `max_value` from the remaining rows. Keep the first three columns
in their original order if you also select columns and want the result to
remain a `MizerScan`. Base `subset()` also preserves this metadata for row
selection. Other data-frame operations can drop it or leave derived attributes
stale: `transform()` drops the scan wrapper, while `dplyr::filter()` can retain
it without refreshing the recorded maxima. For arbitrary data manipulation,
use `as.data.frame(scan)` and make your own plot, or rebuild with `MizerScan()`
and the original metadata before using `plot()` on the scan again.

### Comparing forward and backward sweeps

With `continuation = TRUE` (the default), each value starts from the state
reached at the preceding value. This can reveal **hysteresis** when a model has
coexisting attractors: an increasing and a decreasing sweep can follow
different stable branches at the same parameter values.

Start from models settled near the two ends of the range, with initial states
in the basins of the branches you want to follow. Set the endpoint's scanned
parameter and settle it with `findSteadyState()`, which keeps model parameters
fixed. `tuneSteadyState()` recalibrates reproduction and the resource, so
independently tuned endpoints may no longer describe the same model. Here
`params_low` and `params_high` have the same model parameters apart from fishing
effort:

```r
values <- seq(0, 1.5, 0.05)
up <- scanModel(params_low, scan_values = values,
                set_func = scanEffort(), value_func = getSSB,
                continuation = TRUE, current_scan_value = NULL)
down <- scanModel(params_high, scan_values = rev(values),
                  set_func = scanEffort(), value_func = getSSB,
                  continuation = TRUE, current_scan_value = NULL)
plot(up, style = "envelope")
plot(down, style = "envelope")
```

Leaving `current_scan_value = NULL` preserves the supplied scan order. Giving
it a value, including `"auto"`, instead sweeps outwards from that value,
restarting from the supplied model in each direction. `continuation = FALSE`
restarts every point from the supplied initial
model, so it does not follow the preceding branch.

An abrupt jump between branches can indicate a saddle-node fold. Check
`attractor`, `termination` and `residual` near the jump, allow enough time to
settle, and refine `scan_values` there before interpreting it as a threshold.
These are sweeps of the dynamics: they follow attracting states, not the
unstable branch between them, and do not locate a fold exactly. The `params`
attribute of a scan holds its **starting model**, not the final state of the
sweep; it cannot be used as the settled starting state for a reverse sweep.

## Numerical caveat

`getStability()` involves no time step, so the temporal numerical diffusion of the implicit solver does not enter its answer at all. But the *spatial* numerical diffusion from the default first-order upwind scheme does: it is part of the semi-discretised model that the eigenvalues describe. A real limit cycle can still be damped to a flat line by the spatial scheme alone!

To accurately simulate the fully nonlinear oscillation and confirm a stable cycle, build the model with `second_order_w = TRUE` and project with `method = "second_order"`. See the `run-simulation` skill.
