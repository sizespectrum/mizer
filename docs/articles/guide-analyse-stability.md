# Guide: Analysing dynamic stability

This guide gives an overview of the experimental tools for analysing the
dynamic stability of a mizer steady state and characterising the limit
cycles that can replace it. For full documentation of each function,
follow the links.

Tools for asking whether a mizer steady state is dynamically **stable**
— and, when it is not, characterising the **limit cycle** that replaces
it. These are **experimental**: their interface may still change.
`findSteadyState(solver = "newton")` assumes the standard semichemostat
resource dynamics and solves for the resource alongside the fish, so the
resource density and the feeding levels it implies are self-consistent
even where consumers are satiated.
[`getStability()`](https://sizespectrum.org/mizer/reference/getStability.md)
and
[`getDiscreteStability()`](https://sizespectrum.org/mizer/reference/getDiscreteStability.md)
perturb the resource alongside the fish and so work with any resource
dynamics.

Distinct from calibration: the [guide to reaching steady state and
calibrating](https://sizespectrum.org/mizer/articles/guide-calibrate-model.md)
gets you *onto* a fixed point; what follows analyses the *dynamics
around* it. The natural entry point is
[`findSteadyState(params, solver = "newton")`](https://sizespectrum.org/mizer/reference/findSteadyState.md),
which converges even when the fixed point is dynamically unstable,
unlike the default `solver = "project"`.

------------------------------------------------------------------------

## Is the steady state stable? — `getStability()`

`getStability(params)` computes the continuous-time eigenvalues of the
model. Mizer discretises the size axis but not time, so on the size grid
the model is a system of ODEs;
[`getStability()`](https://sizespectrum.org/mizer/reference/getStability.md)
differentiates their right-hand side directly and takes the eigenvalues
of that Jacobian. **No time step is involved**, so the eigenvalues are a
property of the model rather than of a solver. The result reports
whether the state is stable or unstable, the maximum real part of the
eigenvalues (`max_real_part`; \> 0 means unstable), and the **period**
at which the model rings (`dominant_period` for the dominant mode,
`oscillation_period` for the leading *oscillatory* one). A complex pair
is an oscillatory mode, not a bifurcation: its period is the period of
an emerging limit cycle only where its real part is zero, which is a
Hopf bifurcation and which a single spectrum cannot show. Watch
`Re(leading_oscillatory_eigenvalue)` cross zero over a
[`scanModel()`](https://sizespectrum.org/mizer/reference/scanModel.md)
before calling one. It also returns `leading_eigenvectors`, the top two
eigenvectors of the state space split into `$fish`, a complex array
`(n_species, n_sizes, 2)`, and `$resource`, a complex matrix
`(n_w_full, 2)`, each normalised to maximum modulus 1 (the spatial shape
of the oscillation).

``` r

params <- findSteadyState(params, solver = "newton")  # sit exactly on the (possibly unstable) fixed point
stab   <- getStability(params)
stab                                     # stable/unstable, growth rate, cycle period
```

- `effort` sets the fishing effort used when forming the Jacobian.
- `h` is the relative finite-difference step. Re-running with a
  different `h` is the cheapest check that the model is smooth enough
  for the analysis to mean anything: if the answer moves, do not trust
  it.

**Both
[`getStability()`](https://sizespectrum.org/mizer/reference/getStability.md)
and
[`getOscillationModeSim()`](https://sizespectrum.org/mizer/reference/getOscillationModeSim.md)
linearise at the state stored in the object**, so a model that is not on
a fixed point gets eigenvalues for the neighbourhood of a point it is
not sitting at. Both now warn when handed one; the fix is to run
[`findSteadyState()`](https://sizespectrum.org/mizer/reference/findSteadyState.md)
(or
[`tuneSteadyState()`](https://sizespectrum.org/mizer/reference/tuneSteadyState.md))
first, not to ignore the warning.
[`plot(getSteadyResidual(params))`](https://sizespectrum.org/mizer/reference/plot.md)
shows how far off it is — see the [guide to reaching steady state and
calibrating](https://sizespectrum.org/mizer/articles/guide-calibrate-model.md).

[`tuneSteadyState()`](https://sizespectrum.org/mizer/reference/tuneSteadyState.md)
and
[`findSteadyState()`](https://sizespectrum.org/mizer/reference/findSteadyState.md)
attach a related `"convergence"` attribute, which answers three
questions separately: `termination` says why the run stopped,
`converged` whether the solver met its own criterion, and `attractor`
what the state reached actually is — `"fixed_point"`, `"limit_cycle"` or
`NA`. **Read `attractor`**: it is the only one set from the measured
biomass drift (reported as `residual`), and so the only one that can be
taken as a claim that the model is at a steady state.

Stopping on a fixed point takes both a distance function below
`distance_tol` and a drift within `residual_tol`; a state that meets the
first but not the second keeps running. Limit cycles are detected from
the per-species biomass series sampled at every time step, at every
check rather than only when the distance criterion has failed, because a
cycle whose period divides `t_check` is sampled at one phase and looks
converged. The relative-amplitude floor for calling an oscillation a
cycle is the `amplitude_tol` argument (default `0.01`).

### What the analysis covers

[`getStability()`](https://sizespectrum.org/mizer/reference/getStability.md),
[`getDiscreteStability()`](https://sizespectrum.org/mizer/reference/getDiscreteStability.md),
[`tuneSteadyState()`](https://sizespectrum.org/mizer/reference/tuneSteadyState.md)
and `findSteadyState(solver = "newton")` all treat components registered
with
[`setComponent()`](https://sizespectrum.org/mizer/reference/setComponent.md)
as fixed inputs: they are held at their stored values and have no rows
in the Jacobian. The verdict is therefore about the consumer-resource
subsystem with the components frozen, which is exact when a component is
a fixed input and a good approximation when it is much faster or much
slower than the fish. Mizer warns when it meets a component with
dynamics of its own. `findSteadyState(solver = "project")` is the one
that advances everything.

------------------------------------------------------------------------

## Is mizer’s time step stable? — `getDiscreteStability()`

A different question, and the one to ask when a simulation disagrees
with
[`getStability()`](https://sizespectrum.org/mizer/reference/getStability.md).
`getDiscreteStability(params, dt = 1)` linearises the one-step map that
[`project(method = "euler")`](https://sizespectrum.org/mizer/reference/project.md)
takes at that step size and returns its `discrete_eigenvalues` and
`spectral_radius`. Below 1 means the *numerical scheme* does not amplify
perturbations at that `dt`. It is that step and not an approximation to
it: the same transport coefficients, the same Thomas solver and the
model’s own
[`resource_dynamics`](https://sizespectrum.org/mizer/reference/setResource.md).

The two verdicts can differ, and that is the point: the implicit
transport solve damps oscillations artificially, so an unstable steady
state can have a spectral radius below 1 at a large `dt` and a
simulation will then settle onto a state the model does not hold. If
[`getStability()`](https://sizespectrum.org/mizer/reference/getStability.md)
says unstable but the simulation goes flat, check
[`getDiscreteStability()`](https://sizespectrum.org/mizer/reference/getDiscreteStability.md)
at the `dt` you are projecting with, and reduce `dt` (or use
`method = "tr_bdf2"`).

The discrete eigenvalues cannot be converted into continuous ones by any
exact algebraic relation, because mizer’s step is not fully implicit:
the rates are evaluated at the state at the start of the step and only
the transport solve is implicit. That is why
[`getStability()`](https://sizespectrum.org/mizer/reference/getStability.md)
differentiates the rates of change directly instead of inverting the
one-step map.

------------------------------------------------------------------------

## Visualising the oscillation — `getOscillationModeSim()`

`getOscillationModeSim(params)` takes the output of
[`findSteadyState()`](https://sizespectrum.org/mizer/reference/findSteadyState.md)
and builds a
[`MizerSim`](https://sizespectrum.org/mizer/reference/MizerSim.md)
covering **exactly one period** of the leading oscillatory mode in the
linear approximation,

\[ x(t) = x^\* + A,!, \]

over the whole state (x = (N, n\_{pp})), where (v) is
[`getStability()`](https://sizespectrum.org/mizer/reference/getStability.md)’s
`leading_oscillatory_eigenvector` — the eigenvector of the dominant
*oscillatory* mode, which is not in general the dominant mode — and (A)
is scaled so that the largest relative swing in **species biomass**
equals the `amplitude` argument (default 10%): the hardest-oscillating
species departs that far from its steady biomass and no species departs
further. Fish and resource share one (A), so the resource oscillates
with the phase the mode gives it; on
[`NS_params`](https://sizespectrum.org/mizer/reference/NS_params.md) at
effort 1.5 it leads the fish biomass by about 0.3 years of the 5.12-year
period. The growth factor (e^{t}) is deliberately omitted so the
oscillation closes after one period; the run ends exactly at (T), with a
shortened final interval when `t_save` does not divide it, so the last
state is the first state again. () is recorded as `growth_rate` in the
result’s `sim_params`, and it is what decides whether this is a cycle
the model settles onto or a transient it rings with on the way back.

A biomass amplitude does not bound the individual size classes — a
cohort trough is a far larger relative excursion than the biomass
integral over it — so
[`getOscillationModeSim()`](https://sizespectrum.org/mizer/reference/getOscillationModeSim.md)
reports when the perturbation drives an abundance negative and it has to
be clipped at zero. On `NS_params` at effort 1.5 that starts at around
`amplitude = 0.2`. Past that point the picture is no longer the linear
mode. The result is an ordinary `MizerSim`, so plot it with the standard
tools (see the [guide to analysing and plotting mizer
results](https://sizespectrum.org/mizer/articles/guide-analyse-and-plot.md)):

``` r

params <- findSteadyState(params, solver = "newton")
sim    <- getOscillationModeSim(params, amplitude = 0.1)
plotBiomass(sim)                         # biomass oscillation over one period
animate(plotSpectra(sim))               # the travelling wave in the spectrum
```

This is the *linear* mode — the shape near onset. For the fully
nonlinear cycle, project the real dynamics with the second-order scheme
(see the [guide to running a mizer
simulation](https://sizespectrum.org/mizer/articles/guide-run-simulation.md)).

------------------------------------------------------------------------

## Where does the cycle appear? — `scanModel()`

[`scanModel()`](https://sizespectrum.org/mizer/reference/scanModel.md)
sweeps any aspect of the model and, at each value, follows the attractor
of the full dynamics and records the long-term **range** of a summary
quantity. Draw that range with `style = "envelope"` and you have a
bifurcation diagram: a stable steady state shows as a single line, a
limit cycle as a **band** between its minimum and maximum, so a Hopf
bifurcation appears as the value at which the band opens up.

``` r

scan <- scanModel(params, scan_values = seq(0, 1.5, 0.05),
                  set_func = scanEffort(), value_func = getYield)
plot(scan, style = "envelope")
```

- `set_func` says what to vary —
  [`scanEffort()`](https://sizespectrum.org/mizer/reference/scanEffort.md)
  for the fishing effort of every gear together,
  [`scanFishingMortality()`](https://sizespectrum.org/mizer/reference/scanEffort.md)
  for the fishing mortality on one species,
  [`scanSpeciesParam()`](https://sizespectrum.org/mizer/reference/scanEffort.md)
  for any species parameter. `value_func` says what to measure: any of
  [`getBiomass()`](https://sizespectrum.org/mizer/reference/getBiomass.md),
  [`getYield()`](https://sizespectrum.org/mizer/reference/getYield.md),
  [`getSSB()`](https://sizespectrum.org/mizer/reference/getSSB.md),
  [`getN()`](https://sizespectrum.org/mizer/reference/getN.md).
  `species` restricts which series are shown. See the [guide to
  analysing and plotting mizer
  results](https://sizespectrum.org/mizer/articles/guide-analyse-and-plot.md).
- The settling stage runs
  [`projectUntilSettled()`](https://sizespectrum.org/mizer/reference/projectUntilSettled.md),
  whose `distance_tol`, `residual_tol`, `amplitude_tol` and
  `extinction_threshold` are exposed for tuning when a run is slow to
  settle or a species is collapsing.
- The returned
  [`MizerScan`](https://sizespectrum.org/mizer/reference/MizerScan.md)
  object is a data frame, so the numbers behind the diagram are
  available without re-running the scan.

------------------------------------------------------------------------

## Numerical caveat

[`getStability()`](https://sizespectrum.org/mizer/reference/getStability.md)
involves no time step, so the temporal numerical diffusion of the
implicit solver does not enter its answer at all. But the *spatial*
numerical diffusion from the default first-order upwind scheme does: it
is part of the semi-discretised model that the eigenvalues describe. A
real limit cycle can still be damped to a flat line by the spatial
scheme alone!

To accurately simulate the fully nonlinear oscillation and confirm a
stable cycle, build the model with `second_order_w = TRUE` and project
with `method = "tr_bdf2"`. See the [guide to running a mizer
simulation](https://sizespectrum.org/mizer/articles/guide-run-simulation.md).
