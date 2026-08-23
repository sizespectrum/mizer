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
  found is stable.
---

# Analysing dynamic stability

Tools for asking whether a mizer steady state is dynamically **stable** — and,
when it is not, characterising the **limit cycle** that replaces it. These are
**experimental**: their interface may still change. `findSteadyState(solver =
"newton")` solves for the resource alongside the fish, so the resource density
and the feeding levels it implies are self-consistent even where consumers are
satiated. It works with the built-in semichemostat and logistic dynamics; the
logistic solve currently covers only the positive-resource branch.
`getStability()` and `getDiscreteStability()` perturb the resource alongside the
fish and so work with any resource dynamics.

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
reports whether the state is stable or unstable, the maximum real part of the
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

**Both `getStability()` and `getOscillationModeSim()` linearise at the state stored in
the object**, so a model that is not on a fixed point gets eigenvalues for the
neighbourhood of a point it is not sitting at. Both now warn when handed one;
the fix is to run `findSteadyState()` (or `tuneSteadyState()`) first, not to
ignore the
warning. `plot(getSteadyResidual(params))` shows how far off it is — see the
`calibrate-model` skill.

`tuneSteadyState()` and `findSteadyState()` attach a related `"convergence"`
attribute, which answers three questions separately: `termination` says why the
run stopped, `converged` whether the solver met its own criterion, and
`attractor` what the state reached actually is — `"fixed_point"`,
`"limit_cycle"` or `NA`. **Read `attractor`**: it is the only one set from the
measured biomass drift (reported as `residual`), and so the only one that can be
taken as a claim that the model is at a steady state.

Stopping on a fixed point takes both a distance function below `distance_tol`
and a drift within `residual_tol`; a state that meets the first but not the
second keeps running. Limit cycles are detected from the per-species biomass
series sampled at every time step, at every check rather than only when the
distance criterion has failed, because a cycle whose period divides `t_check` is
sampled at one phase and looks converged. The relative-amplitude floor for calling an
oscillation a cycle is the `amplitude_tol` argument (default `0.01`).

### What the analysis covers

`getStability()`, `getDiscreteStability()`, `tuneSteadyState()` and
`findSteadyState(solver = "newton")` all treat components registered with
`setComponent()` as fixed inputs: they are held at their stored values and have
no rows in the Jacobian. The verdict is therefore about the consumer-resource
subsystem with the components frozen, which is exact when a component is a fixed
input and a good approximation when it is much faster or much slower than the
fish. Mizer warns when it meets a component with dynamics of its own.
`findSteadyState(solver = "project")` is the one that advances everything.

## Is mizer's time step stable? — `getDiscreteStability()`

A different question, and the one to ask when a simulation disagrees with
`getStability()`. `getDiscreteStability(params, dt = 1)` linearises the one-step
map that `project(method = "euler")` takes at that step size and returns its
`discrete_eigenvalues` and `spectral_radius`. Below 1 means the *numerical
scheme* does not amplify perturbations at that `dt`. It is that step and not an
approximation to it: the same transport coefficients, the same Thomas solver and
the model's own `resource_dynamics`.

The two verdicts can differ, and that is the point: the implicit transport solve
damps oscillations artificially, so an unstable steady state can have a spectral
radius below 1 at a large `dt` and a simulation will then settle onto a state the
model does not hold. If `getStability()` says unstable but the simulation goes
flat, check `getDiscreteStability()` at the `dt` you are projecting with, and
reduce `dt` (or use `method = "tr_bdf2"`).

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

## Numerical caveat

`getStability()` involves no time step, so the temporal numerical diffusion of the implicit solver does not enter its answer at all. But the *spatial* numerical diffusion from the default first-order upwind scheme does: it is part of the semi-discretised model that the eigenvalues describe. A real limit cycle can still be damped to a flat line by the spatial scheme alone!

To accurately simulate the fully nonlinear oscillation and confirm a stable cycle, build the model with `second_order_w = TRUE` and project with `method = "tr_bdf2"`. See the `run-simulation` skill.
