# Superseded names for the steady-state finders

**\[superseded\]**

Neither of these names says what distinguishes the two functions, and
`projectToSteady()` returns a different class depending on an argument.
Both have been replaced:

|                     |                                                                                                                                                                                                      |
|---------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Superseded          | Use instead                                                                                                                                                                                          |
| `steady()`          | [`tuneSteadyState()`](https://sizespectrum.org/mizer/reference/tuneSteadyState.md)                                                                                                                   |
| `projectToSteady()` | [`findSteadyState()`](https://sizespectrum.org/mizer/reference/findSteadyState.md), or [`projectUntilSettled()`](https://sizespectrum.org/mizer/reference/projectUntilSettled.md) for the trajectory |

The distinction the new names carry is what each one holds fixed.
[`tuneSteadyState()`](https://sizespectrum.org/mizer/reference/tuneSteadyState.md)
holds the inputs to the fish dynamics — the reproduction rate and the
resource — at the values you supply while the spectra settle, and then
adjusts the parameters that generate them (`erepro`/`R_max` and `cc_pp`)
so that those held values are steady too. That is what `steady()` always
did.
[`findSteadyState()`](https://sizespectrum.org/mizer/reference/findSteadyState.md)
changes no parameter and lets reproduction, the resource and the spectra
settle together, which is what `projectToSteady()` did.

Each of the two also gained a `solver` argument, so the same job can be
done either by running the dynamics (`solver = "project"`, the default
and the old behaviour) or by solving the steady-state equation directly
with a Newton-type root finder (`solver = "newton"`), which converges
even at a dynamically unstable steady state.

The `return_sim` argument is gone from the new functions:
[`tuneSteadyState()`](https://sizespectrum.org/mizer/reference/tuneSteadyState.md)
and
[`findSteadyState()`](https://sizespectrum.org/mizer/reference/findSteadyState.md)
always return a `MizerParams` object and
[`projectUntilSettled()`](https://sizespectrum.org/mizer/reference/projectUntilSettled.md)
always returns a `MizerSim`.

The old names are however not going away. They are thin wrappers that
reproduce the old behaviour exactly, `return_sim` included: they do not
warn and they will keep working, so existing code and old scripts run
unchanged. They are not used anywhere inside mizer and should not be
used in new code.

## Usage

``` r
steady(
  params,
  t_max = 100,
  t_per = 1.5,
  dt = 0.1,
  t_save = dt,
  tol = 0.1 * dt,
  amplitude_tol = 0.01,
  amp_rel_tol = 0.01,
  extinction_threshold = 1e-06,
  return_sim = FALSE,
  preserve = c("reproduction_level", "erepro", "R_max"),
  progress_bar = TRUE,
  info_level = default_info_level(),
  method = c("euler", "predictor_corrector", "tr_bdf2")
)

projectToSteady(
  params,
  effort = params@initial_effort,
  distance_func = distanceSSLogN,
  t_per = 1.5,
  t_max = 100,
  dt = 0.1,
  t_save = dt,
  tol = 0.1 * t_per,
  amplitude_tol = 0.01,
  amp_rel_tol = 0.1,
  extinction_threshold = 1e-06,
  return_sim = FALSE,
  progress_bar = TRUE,
  info_level = default_info_level(),
  method = c("euler", "predictor_corrector", "tr_bdf2"),
  ...
)
```

## Arguments

- params:

  A
  [MizerParams](https://sizespectrum.org/mizer/reference/MizerParams-class.md)
  object

- t_max:

  The maximum number of years to run the simulation. Default is 100.

- t_per:

  The interval in years at which convergence is checked, and hence also
  the interval at which the trajectory is saved when
  `return_sim = TRUE`. In
  [`projectUntilSettled()`](https://sizespectrum.org/mizer/reference/projectUntilSettled.md)
  these two roles have been separated into `t_check` and `t_save`.

- dt:

  The time step to use in
  [`project()`](https://sizespectrum.org/mizer/reference/project.md).

- t_save:

  Has no effect. It briefly controlled how finely the biomass series
  used for limit-cycle detection was sampled; that series is now sampled
  at every time step, which is what its default `dt` gave.

- tol:

  The simulation stops when the relative change in the egg production
  RDI over t_per years is less than tol for every species.

- amplitude_tol:

  **\[experimental\]** The minimum relative biomass amplitude for a
  persistent oscillation to be reported as a limit cycle rather than
  treated as an (effectively steady) fixed point. This is a fraction of
  mean biomass and is kept separate from `distance_tol` (which measures
  convergence to a fixed point on a different scale). Default `0.01`.

- amp_rel_tol:

  **\[experimental\]** Maximum relative change of amplitude between
  successive periods for the cycle to count as settled. Default `0.01`.

- extinction_threshold:

  **\[experimental\]** A species is treated as going extinct, stopping
  the run, once its reproduction rate (RDD) falls below this fraction of
  its value at the start of the run. For example the default `1e-6`
  treats a species as extinct once its reproduction has collapsed to a
  millionth of its initial level. Because it is relative to the initial
  reproduction, a species that starts with zero reproduction is flagged
  immediately, and (in
  [`tuneSteadyState()`](https://sizespectrum.org/mizer/reference/tuneSteadyState.md),
  where reproduction is held constant) a healthy species is never
  flagged.

- return_sim:

  If TRUE, the function returns the MizerSim object holding the result
  of the simulation run, saved at intervals of `t_per`. If FALSE
  (default) the function returns a MizerParams object with the "initial"
  slots set to the steady state.

- preserve:

  **\[experimental\]** Specifies whether the `reproduction_level` should
  be preserved (default) or the maximum reproduction rate `R_max` or the
  reproductive efficiency `erepro`. See
  [`setBevertonHolt()`](https://sizespectrum.org/mizer/reference/setBevertonHolt.md)
  for an explanation of the `reproduction_level`.

- progress_bar:

  A shiny progress object to implement a progress bar in a shiny app.
  Default FALSE.

- info_level:

  Controls the amount of information messages that are shown. Higher
  levels lead to more messages, `info_level = 0` gives silence. The
  default is taken from the `mizer_info_level` option, see
  [`default_info_level()`](https://sizespectrum.org/mizer/reference/default_info_level.md).

- method:

  The numerical method to use for the consumer density update. See
  [`project()`](https://sizespectrum.org/mizer/reference/project.md).

- effort:

  The fishing effort to use throughout. By default the initial effort
  stored in `params`.

- distance_func:

  A function that will be called at every check with both the previous
  and the new state and that should return a number that in some sense
  measures the distance between the states. By default this uses the
  function
  [`distanceSSLogN()`](https://sizespectrum.org/mizer/reference/distanceSSLogN.md)
  that you can use as a model for your own distance function.

- ...:

  Further arguments will be passed on to your distance function.

## Value

A `MizerParams` object, or a `MizerSim` object if `return_sim = TRUE`,
in either case carrying the `"convergence"` attribute described in
[`projectUntilSettled()`](https://sizespectrum.org/mizer/reference/projectUntilSettled.md).

## See also

[`tuneSteadyState()`](https://sizespectrum.org/mizer/reference/tuneSteadyState.md),
[`findSteadyState()`](https://sizespectrum.org/mizer/reference/findSteadyState.md),
[`projectUntilSettled()`](https://sizespectrum.org/mizer/reference/projectUntilSettled.md)
