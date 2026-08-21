---
name: analyse-stability
description: >-
  Analyse the dynamic stability of a mizer steady state and characterise limit
  cycles (experimental). Use whenever the user asks whether a steady state is
  stable or unstable, wants the leading eigenvalue, the period of an emergent
  oscillation, a Hopf bifurcation, the spectral radius of the numerical time
  step, a limit cycle to build or plot, or a bifurcation diagram over fishing
  effort — via getStability(), getDiscreteStability(), steadyNewton(),
  getLimitCycleSim() and scanModel(). This
  skill and calibrate-model share steadyNewton() and getSteadyResidual(): use
  calibrate-model to find a steady state, this skill to ask whether the state you
  found is stable. Assumes the standard semichemostat resource dynamics.
---

# Analysing dynamic stability

Tools for asking whether a mizer steady state is dynamically **stable** — and,
when it is not, characterising the **limit cycle** that replaces it. These are
**experimental**: their interface may still change. All assume the standard
semichemostat resource dynamics. `steadyNewton()` solves for the resource
alongside the fish, so the resource density and the feeding levels it implies
are self-consistent even where consumers are satiated; `getStability()` instead
treats the resource as a quasi-static fast variable, and
`include_resource = TRUE` turns that approximation off.

Distinct from calibration: the `calibrate-model` skill gets you *onto* a fixed
point; what follows analyses the *dynamics around* it. The natural entry point
is `steadyNewton()`, which converges even when the fixed point is dynamically
unstable, unlike `steady()`.

## Is the steady state stable? — `getStability()`

`getStability(params)` computes the continuous-time eigenvalues of the model.
Mizer discretises the size axis but not time, so on the size grid the model is a
system of ODEs; `getStability()` differentiates their right-hand side directly
and takes the eigenvalues of that Jacobian. **No time step is involved**, so the
eigenvalues are a property of the model rather than of a solver. The result
reports whether the state is stable or unstable, the maximum real part of the
eigenvalues (`max_real_part`; > 0 means unstable), and — when the system is near
a Hopf bifurcation — the **period** of the emergent limit cycle
(`dominant_period`, `hopf_period`). It also returns `leading_eigenvectors`, a
complex array `(n_species, n_sizes, 2)` of the top two eigenvectors in
fish-abundance space, normalised to maximum modulus 1 (the spatial shape of the
oscillation).

```r
params <- steadyNewton(params)          # sit exactly on the (possibly unstable) fixed point
stab   <- getStability(params)
stab                                     # stable/unstable, growth rate, cycle period
```

- `include_resource = TRUE` computes the full coupled (fish + resource) Jacobian
  instead of the quasi-static approximation — mainly to verify that the
  approximation makes little difference.
- `effort` / `reproduction` set the fishing effort and reproduction handling used
  when forming the Jacobian.
- `h` is the relative finite-difference step. Re-running with a different `h` is
  the cheapest check that the model is smooth enough for the analysis to mean
  anything: if the answer moves, do not trust it.

**Both `getStability()` and `getLimitCycleSim()` linearise at the state stored in
the object**, so a model that is not on a fixed point gets eigenvalues for the
neighbourhood of a point it is not sitting at. Both now warn when handed one;
the fix is to run `steadyNewton()` (or `steady()`) first, not to ignore the
warning. `plot(getSteadyResidual(params))` shows how far off it is — see the
`calibrate-model` skill.

`steady()` and `projectToSteady()` attach a related `"convergence"` attribute
recording whether the run **dropped below its tolerance** (`type =
"below_tolerance"`, which suggests but does not prove a fixed point), settled on
**a limit cycle**, or **neither**, together with the cycle period and relative
amplitude when a cycle is found. Limit cycles are detected from the per-species biomass series sampled at
`t_save`; the relative-amplitude floor for calling an oscillation a cycle is the
`amplitude_tol` argument (default `0.01`), independent of the fixed-point
tolerance `tol`.

## Is mizer's time step stable? — `getDiscreteStability()`

A different question, and the one to ask when a simulation disagrees with
`getStability()`. `getDiscreteStability(params, dt = 1)` linearises the one-step
map that `project(method = "euler")` takes at that step size and returns its
`discrete_eigenvalues` and `spectral_radius`. Below 1 means the *numerical
scheme* does not amplify perturbations at that `dt`.

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

## Visualising the limit cycle — `getLimitCycleSim()`

`getLimitCycleSim(params)` takes the output of `steadyNewton()` and builds a
`MizerSim` covering **one period** of the limit cycle in the linear approximation,

\[ N(t) = N^* + A\,\mathrm{Re}\!\left[e^{i\theta t}\,\mathbf v\right], \]

where \(\mathbf v\) is the leading complex eigenvector and the amplitude \(A\) is
scaled so the maximum relative perturbation equals the `amplitude` argument
(default 10%). The result is an ordinary `MizerSim`, so plot it with the standard
tools (see the `analyse-and-plot` skill):

```r
params <- steadyNewton(params)
sim    <- getLimitCycleSim(params, amplitude = 0.1)
plotBiomass(sim)                         # biomass oscillation over one period
animate(plotSpectra(sim))               # the travelling wave in the spectrum
```

This is the *linear* cycle — the shape near onset. For the fully nonlinear cycle,
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
- The settling stage runs `projectToSteady()`, whose `tol`, `amplitude_tol` and
  `extinction_threshold` are exposed for tuning when a run is slow to settle or a
  species is collapsing.
- The returned `MizerScan` object is a data frame, so the numbers behind the
  diagram are available without re-running the scan.

## Numerical caveat

`getStability()` involves no time step, so the temporal numerical diffusion of the implicit solver does not enter its answer at all. But the *spatial* numerical diffusion from the default first-order upwind scheme does: it is part of the semi-discretised model that the eigenvalues describe. A real limit cycle can still be damped to a flat line by the spatial scheme alone!

To accurately simulate the fully nonlinear oscillation and confirm a stable cycle, build the model with `second_order_w = TRUE` and project with `method = "tr_bdf2"`. See the `run-simulation` skill.
