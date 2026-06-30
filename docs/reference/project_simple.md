# Project abundances by a given number of time steps into the future

This is an internal function used by the user-facing
[`project()`](https://sizespectrum.org/mizer/reference/project.md)
function. It is of potential interest only to mizer extension authors.

## Usage

``` r
project_simple(
  params,
  n,
  n_pp,
  n_other,
  effort,
  t,
  dt,
  steps,
  resource_dynamics_fn,
  other_dynamics_fns,
  rates_fns,
  method = c("euler", "predictor_corrector", "tr_bdf2"),
  ...
)
```

## Arguments

- params:

  A MizerParams object.

- n:

  An array (species x size) with the number density at start of
  simulation.

- n_pp:

  A vector (size) with the resource number density at start of
  simulation.

- n_other:

  A named list with the abundances of other components at start of
  simulation.

- effort:

  The fishing effort to be used throughout the simulation. This must be
  a vector or list with one named entry per fishing gear.

- t:

  Time at the start of the simulation.

- dt:

  Size of time step.

- steps:

  The number of time steps by which to project.

- resource_dynamics_fn:

  The function for the resource dynamics. See Details.

- other_dynamics_fns:

  List with the functions for the dynamics of the other components. See
  Details.

- rates_fns:

  List with the functions for calculating the rates. See Details.

- method:

  The numerical method to use for the consumer density update. See
  [`project()`](https://sizespectrum.org/mizer/reference/project.md).

- ...:

  Other arguments that are passed on to the rate functions.

## Value

List with the final values of `n`, `n_pp`, and `n_other`, together with
`rates`, the rates calculated at the start of the final update step.

## Details

The function does not check its arguments because it is meant to be as
fast as possible to allow it to be used in a loop. For example, it is
called in
[`project()`](https://sizespectrum.org/mizer/reference/project.md) once
for every saved value. The function also does not save its intermediate
results but only returns the result at time `t + dt * steps`. During
this time it uses the constant fishing effort `effort`.

The functional arguments can be calculated from slots in the `params`
object with

    resource_dynamics_fn <- get(params@resource_dynamics)
    other_dynamics_fns <- lapply(params@other_dynamics, get)
    rates_fns <- lapply(params@rates_funcs, get)

The reason the function does not do that itself is to shave 20
microseconds of its running time, which pays when the function is called
hundreds of times in a row.

This function is also used in
[`steady()`](https://sizespectrum.org/mizer/reference/steady.md). In
between calls to `project_simple()` the
[`steady()`](https://sizespectrum.org/mizer/reference/steady.md)
function checks whether the values are still changing significantly, so
that it can stop when a steady state has been approached. Mizer
extension packages might have a similar need to run a simulation
repeatedly for short periods to run some other code in between. Because
this code may want to use the values of the rates from the final update
step, these too are included in the returned list.
