---
name: run-simulation
description: >-
  Project a mizer model forward in time and set up fishing-effort scenarios. Use
  whenever the user wants to run a simulation with project(), specify constant or
  time-varying fishing effort, choose a projection method or time step, run to a
  new steady state after a change, continue an existing MizerSim, or set up
  scenario comparisons. For extracting and plotting the results, see the
  analyse-and-plot skill.
---

# Running a mizer simulation

`project()` advances a `MizerParams` object through time and returns a
`MizerSim`. The params object must already be set up and (usually) at steady
state — see the `build-multispecies-model` skill and the `calibrate-model` skill.

```r
sim <- project(params, t_max = 20, effort = 1)
```

## Key `project()` arguments

| Argument | Meaning |
|---|---|
| `object` | a `MizerParams` (fresh run) or a `MizerSim` (continue from its end) |
| `effort` | fishing effort — see the four forms below |
| `t_max` | number of years to simulate (default 100) |
| `dt` | integration time step (default 0.1); reduce if the run is unstable |
| `t_save` | interval (years) at which output is stored (default 1) |
| `t_start` | initial time / calendar year for the output (default 0) |
| `method` | `"euler"` (default), `"predictor_corrector"`, or `"tr_bdf2"` |
| `callback` | a function called at each saved step (e.g. to log or intervene) |
| `progress_bar` | set `FALSE` to silence the progress bar |

## Specifying fishing effort

`effort` can be given four ways. If omitted, the model's stored
`initial_effort` is used.

```r
project(params, effort = 1)                          # 1. scalar: same for all gears, constant
project(params, effort = c(Otter = 0.5, Beam = 1))   # 2. named vector: per-gear, constant
project(params, effort = c(0.5, 1, 0))               # 3. vector in gear order, constant
project(params, effort = effort_array)               # 4. time × gear array: effort through time
```

For form 4, build a `time × gear` matrix with **numeric, increasing** row names
(times) and column names matching the gear names:

```r
gears <- names(initial_effort(params))               # gear names (a named vector)
years <- 2010:2030
effort_array <- array(1, dim = c(length(years), length(gears)),
                      dimnames = list(time = years, gear = gears))
effort_array[as.character(2020:2030), "Otter"] <- 1.5  # ramp one gear up from 2020
sim <- project(params, effort = effort_array)
```

Each effort value applies from its time until the next time in the array. With
an array, the simulation starts at the smallest time; use `t_max` to extend
beyond the last row.

To change *which* gears exist, or their selectivity and catchability, before
running, see the `set-up-fishing` skill.

## Common patterns

**Run to a new steady state after a change.** To get the equilibrium a change
implies, rather than a fixed number of years, use `steady()` or
`projectToSteady()` from the `calibrate-model` skill instead of a long
`project()`.

**Continue a simulation.** Pass a `MizerSim` back to `project()`; it resumes
from the final state:

```r
sim2 <- project(sim, t_max = 10, effort = 2)
```

**Carry a simulation's state into a fresh `MizerParams`** — to start a new run
under different settings, or to analyse a particular time step:

```r
params_end <- finalParams(sim)      # state at the last time step
params_0   <- initialParams(sim)    # state at the first time step
params_t   <- getParams(sim, time_range = 2010:2015)  # averaged over a range
sim_next   <- project(params_end, t_max = 20, effort = 0.5)
```

`finalParams(sim)` and `initialParams(sim)` are shorthands for the final and
initial steps of `getParams()`. Prefer these over the **deprecated**
`setInitialValues(params, sim)`.

**Scenario comparison.** Project the same params under different efforts, then
compare with the plotting tools (see the `analyse-and-plot` skill):

```r
sim_low  <- project(params, t_max = 30, effort = 0.5)
sim_high <- project(params, t_max = 30, effort = 1.5)
plotSpectra2(sim_low, sim_high, "F = 0.5", "F = 1.5")
```

## Numerical scheme: watch for numerical diffusion

For steady states and slow biomass or yield trends the defaults are fine. But
the default flux scheme (first-order **upwind**) carries substantial *numerical
diffusion*: it smears out cohorts and travelling waves and can **silently damp a
real oscillation or limit cycle down to a flat line** — with no error to warn
you. This is a correctness issue, not just an accuracy one, so for **any study
of dynamics** (oscillations, cohort or wave structure, generation cycles, real
growth-diffusion effects) switch to the second-order scheme:

- **Space:** build the model with `second_order_w = TRUE` (the van Leer,
  flux-limited scheme) — set it in `newMultispeciesParams(..., second_order_w = TRUE)`,
  or on an existing model with `second_order_w(params) <- TRUE`. Because it
  changes the discrete steady state it lives in the params object, not in a
  `project()` argument.
- **Time:** project with `method = "tr_bdf2"` (L-stable, second order in time).

```r
params <- newMultispeciesParams(sp, second_order_w = TRUE)
sim    <- project(params, t_max = 200, method = "tr_bdf2")
```

Symptom that you needed this: an oscillation that is present with
`second_order_w = TRUE` but disappears (settles to a flat steady state) under the
default upwind scheme is being killed by numerical diffusion, not by real
dynamics.

**Expect the numbers to move, and don't read that as a bug.** A bare
`second_order_w = TRUE` sets two things: the flux scheme above *and*
`bin_average`, which switches every integral over the size grid from a left-edge
Riemann sum to a proper bin average. Unnormalised quantities — biomass, yield,
`getDiet(proportion = FALSE)` — therefore shift by roughly `(1 + beta) / 2`,
about 10% for a typical grid, where `beta` is the ratio between neighbouring
grid points. That is the more accurate integral, not a regression. Proportions
and ratios are largely unaffected because the factor cancels. Set the two
independently if you want the better flux scheme without moving your summary
numbers:

```r
second_order_w(params) <- c(flux = "van_leer")   # leaves bin_average alone
second_order_w(params)                           # inspect: flux and bin_average
```

The scheme lives in the `MizerParams`, so a `MizerSim` carries it: comparing a
run made under one setting with a run made under the other compares two
discretisations as well as two scenarios. Recalibrate after switching it on.

**Isolating a feedback loop.** To switch off the resource → growth feedback (the
"phantom jam") while keeping everything else — for example to separate an
internal instability from a boundary or reproduction one — freeze the resource
before projecting:

```r
params <- setResource(params, resource_dynamics = "resource_constant")
```

## Tips

- If a run blows up or oscillates unphysically, first reduce `dt` (e.g. `0.01`);
  a stiff model may also do better with `method = "tr_bdf2"`.
- `t_save` controls output resolution, not accuracy — the model always steps at
  `dt` internally.
- The `MizerSim` stores the `MizerParams` it used, so a sim is self-contained
  for later analysis.
