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
  subcycle = FALSE,
  max_courant = 1,
  max_substeps = 10000,
  info_level = default_info_level(),
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

- subcycle:

  Whether to split requested steps using the growth Courant target, as
  in [`project()`](https://sizespectrum.org/mizer/reference/project.md).
  Default `FALSE`.

- max_courant:

  Positive growth Courant target, as in
  [`project()`](https://sizespectrum.org/mizer/reference/project.md).

- max_substeps:

  Positive integer cap per requested step, as in
  [`project()`](https://sizespectrum.org/mizer/reference/project.md).

- info_level:

  Level of information reported; see
  [`default_info_level()`](https://sizespectrum.org/mizer/reference/default_info_level.md).

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

The default method prepares the functional arguments from the model,
including the derivatives of components registered with `rate_fun`.
Callers that repeatedly advance the same model can supply prepared
functions to reuse their dispatch. Component functions must retain the
derivative metadata added during preparation; looking up
`other_dynamics` names alone does not prepare an equation-only component
for second-order integration.

This function is also used by the steady-state finders. In between calls
to `project_simple()` the
[`steady()`](https://sizespectrum.org/mizer/reference/superseded_steady.md)
function checks whether the values are still changing significantly, so
that it can stop when a steady state has been approached. Mizer
extension packages might have a similar need to run a simulation
repeatedly for short periods to run some other code in between. Because
this code may want to use the values of the rates from the final update
step, these too are included in the returned list.
