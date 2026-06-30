# Extract the model state from a simulation

A `MizerParams` object describes the state of the ecosystem: its species
parameters, size grid, rate functions, *and* the current abundances
stored in the `initial_n`, `initial_n_pp`, `initial_n_other`, and
`initial_effort` slots. `getParams()` extracts that state from a
`MizerSim` object, averaged over a chosen time range (or at a single
time point).

## Usage

``` r
getParams(sim, time_range, geometric_mean = FALSE)
```

## Arguments

- sim:

  A `MizerSim` object.

- time_range:

  The time range to average the abundances over. Can be a vector of
  values, a vector of min and max time, or a single value. Only the
  range of times is relevant, i.e., all times between the smallest and
  largest will be selected. Default is the final time step.

- geometric_mean:

  **\[experimental\]** If `TRUE`, the average of the abundances over the
  time range is a geometric mean instead of the default arithmetic mean.
  This does not affect the average of the effort or of other components,
  which is always arithmetic.

## Value

A `MizerParams` object with `initial_n`, `initial_n_pp`,
`initial_n_other`, and `initial_effort` set to the (averaged) values
from the simulation.

## Details

When no `time_range` is given, the state at the final time step is
returned. Use
[`initialParams()`](https://sizespectrum.org/mizer/reference/initialParams.md)
or
[`finalParams()`](https://sizespectrum.org/mizer/reference/finalParams.md)
as convenient shorthand for the state at the initial and final time
respectively.

The abundances set in the returned `MizerParams` object are averages
over the selected time range. By default this is an arithmetic mean; set
`geometric_mean = TRUE` to use a geometric mean instead (this does not
affect the effort or other components, which are always averaged
arithmetically).

## See also

[`initialParams()`](https://sizespectrum.org/mizer/reference/initialParams.md),
[`finalParams()`](https://sizespectrum.org/mizer/reference/finalParams.md)

## Examples

``` r
sim <- project(NS_params, t_max = 20, effort = 0.5)
# Extract state at a specific time
params_2010 <- getParams(sim, time_range = 10)
# Extract state averaged over the last 10 years
params_avg <- getParams(sim, time_range = c(10, 20))
```
