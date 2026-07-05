# Project to steady state

**\[experimental\]**

Run the full dynamics, as in
[`project()`](https://sizespectrum.org/mizer/reference/project.md), but
stop once the change has slowed down sufficiently, in the sense that the
distance between states at successive time steps is less than `tol`. You
determine how the distance is calculated.

## Usage

``` r
projectToSteady(
  params,
  effort = params@initial_effort,
  distance_func = distanceSSLogN,
  t_per = 1.5,
  t_max = 100,
  dt = 0.1,
  tol = 0.1 * t_per,
  return_sim = FALSE,
  progress_bar = TRUE,
  info_level = 3,
  method = c("euler", "predictor_corrector", "tr_bdf2"),
  ...
)
```

## Arguments

- params:

  A
  [MizerParams](https://sizespectrum.org/mizer/reference/MizerParams-class.md)
  object

- effort:

  The fishing effort to be used throughout the simulation. This is
  validated by
  [`validEffortVector()`](https://sizespectrum.org/mizer/reference/validEffortVector.md)
  and can therefore be `NULL`, a single numeric value used for all
  gears, an unnamed numeric vector with one entry per gear, or a named
  numeric vector for some or all gears.

- distance_func:

  A function that will be called after every `t_per` years with both the
  previous and the new state and that should return a number that in
  some sense measures the distance between the states. By default this
  uses the function
  [`distanceSSLogN()`](https://sizespectrum.org/mizer/reference/distanceSSLogN.md)
  that you can use as a model for your own distance function.

- t_per:

  The simulation is broken up into shorter runs of `t_per` years, after
  each of which we check for convergence. Default value is 1.5. This
  should be chosen as an odd multiple of the timestep `dt` in order to
  be able to detect period 2 cycles.

- t_max:

  The maximum number of years to run the simulation. Default is 100.

- dt:

  The time step to use in
  [`project()`](https://sizespectrum.org/mizer/reference/project.md).

- tol:

  The simulation stops when the relative change in the egg production
  RDI over `t_per` years is less than `tol` for every species.

- return_sim:

  If TRUE, the function returns the MizerSim object holding the result
  of the simulation run, saved at intervals of `t_per`. If FALSE
  (default) the function returns a MizerParams object with the "initial"
  slots set to the steady state.

- progress_bar:

  A shiny progress object to implement a progress bar in a shiny app.
  Default FALSE.

- info_level:

  Controls the amount of information messages that are shown. Higher
  levels lead to more messages.

- method:

  The numerical method to use for the consumer density update. See
  [`project()`](https://sizespectrum.org/mizer/reference/project.md).

- ...:

  Further arguments will be passed on to your distance function.

## Value

If `return_sim = FALSE`, a `MizerParams` object with the initial state
replaced by the final state found by the steady-state search. If
`return_sim = TRUE`, a `MizerSim` object containing the intermediate
states saved every `t_per` years.

## See also

[`distanceSSLogN()`](https://sizespectrum.org/mizer/reference/distanceSSLogN.md),
[`distanceMaxRelRDI()`](https://sizespectrum.org/mizer/reference/distanceMaxRelRDI.md)
