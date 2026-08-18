---
name: analyse-stability
description: >-
  Analyse the dynamic stability of a mizer steady state and characterise limit
  cycles (experimental). Use whenever the user asks whether a steady state is
  stable or unstable, wants the spectral radius or leading eigenvalue, the period
  of an emergent oscillation, a Hopf bifurcation, a limit cycle to build or plot,
  or a bifurcation diagram over fishing effort — via getStability(),
  steadyNewton(), getLimitCycleSim() and plotBifurcation(). This
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

`getStability(params)` computes the exact continuous-time eigenvalues (by linearising the one-step map and analytically removing the temporal numerical diffusion). The result reports whether the state is stable or unstable, the maximum real part of the eigenvalues (`max_real_part`; > 0 means unstable), and — when the system is near a Hopf bifurcation — the **period** of the emergent limit cycle. You can optionally pass a `dt` argument to retrieve the `discrete_eigenvalues` and `spectral_radius` of the numerical solver for that exact time step. It also returns `leading_eigenvectors`, a complex array
`(n_species, n_sizes, 2)` of the top two eigenvectors in fish-abundance space,
normalised to maximum modulus 1 (the spatial shape of the oscillation).

```r
params <- steadyNewton(params)          # sit exactly on the (possibly unstable) fixed point
stab   <- getStability(params)
stab                                     # stable/unstable, spectral radius, cycle period
```

- `include_resource = TRUE` computes the full coupled (fish + resource) Jacobian
  instead of the quasi-static approximation — mainly to verify that the
  approximation makes little difference.
- `effort` / `reproduction` set the fishing effort and reproduction handling used
  when forming the map.

**Both `getStability()` and `getLimitCycleSim()` linearise at the state stored in
the object**, so a model that is not on a fixed point gets eigenvalues for the
neighbourhood of a point it is not sitting at. Both now warn when handed one;
the fix is to run `steadyNewton()` (or `steady()`) first, not to ignore the
warning. `plot(getSteadyResidual(params))` shows how far off it is — see the
`calibrate-model` skill.

`steady()` and `projectToSteady()` attach a related `"convergence"` attribute
recording whether the run settled on a **steady state, a limit cycle, or
neither**, together with the cycle period and relative amplitude when a cycle is
found. Limit cycles are detected from the per-species biomass series sampled at
`t_save`; the relative-amplitude floor for calling an oscillation a cycle is the
`amplitude_tol` argument (default `0.01`), independent of the fixed-point
tolerance `tol`.

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

## Where does the cycle appear? — `plotBifurcation()`

`plotBifurcation(params, effort = ...)` draws a bifurcation diagram over fishing
effort. For each effort value it follows the attractor of the full dynamics and
plots the long-term **range** of a summary quantity: a stable steady state shows
as a single line, a limit cycle as a **band** between its minimum and maximum, so
a Hopf bifurcation appears as the effort at which the band opens up.

```r
plotBifurcation(params, effort = seq(0, 1.5, 0.05), value = "yield")
```

- `value` chooses the summary quantity (`"biomass"`, `"yield"`, or `"ssb"`);
  `species` restricts which species are shown.
- The settling stage runs `projectToSteady()`, whose `tol`, `amplitude_tol` and
  `extinction_threshold` are exposed for tuning when a run is slow to settle or a
  species is collapsing.
- `return_data = TRUE` returns the underlying data instead of the plot.

## Numerical caveat

While `getStability()` is now exact for the continuous-time dynamics (by analytically removing the temporal numerical diffusion of the underlying implicit map), the *spatial* numerical diffusion from the default first-order upwind scheme remains. A real limit cycle can still be damped to a flat line by the spatial scheme alone!

To accurately simulate the fully nonlinear oscillation and confirm a stable cycle, build the model with `second_order_w = TRUE` and project with `method = "tr_bdf2"`. See the `run-simulation` skill.
