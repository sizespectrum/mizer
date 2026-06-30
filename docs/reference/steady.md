# Set initial values to a steady state for the model

The steady state is found by running the dynamics while keeping
reproduction, resource and other components constant until the size
spectra no longer change much (or until time `t_max` is reached, if
earlier).

## Usage

``` r
steady(
  params,
  t_max = 100,
  t_per = 1.5,
  dt = 0.1,
  tol = 0.1 * dt,
  return_sim = FALSE,
  preserve = c("reproduction_level", "erepro", "R_max"),
  progress_bar = TRUE,
  info_level = 3,
  method = c("euler", "predictor_corrector", "tr_bdf2")
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

  The simulation is broken up into shorter runs of `t_per` years, after
  each of which we check for convergence. Default value is 1.5. This
  should be chosen as an odd multiple of the timestep `dt` in order to
  be able to detect period 2 cycles.

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
  levels lead to more messages.

- method:

  The numerical method to use for the consumer density update. See
  [`project()`](https://sizespectrum.org/mizer/reference/project.md).

## Value

If `return_sim = FALSE`, a `MizerParams` object with the initial state
replaced by the steady state. If `return_sim = TRUE`, a `MizerSim`
object containing the intermediate states saved every `t_per` years.

## Details

If the model use Beverton-Holt reproduction then the reproduction
parameters are set to values that give the level of reproduction
observed in that steady state. The `preserve` argument can be used to
specify which of the reproduction parameters should be preserved.

## Examples

``` r
# \donttest{
params <- newTraitParams()
species_params(params)$gamma[5] <- 3000
params <- steady(params)
#> Convergence was achieved in 12 years.
plotSpectra(params)

# }
```
