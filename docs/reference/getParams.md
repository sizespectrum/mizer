# Extract the model state from a simulation

A `MizerParams` object describes the state of the ecosystem: its species
parameters, size grid, rate functions, *and* the current abundances
stored in the `initial_n`, `initial_n_pp`, `initial_n_other`, and
`initial_effort` slots. These functions extract that state from a
`MizerSim` object.

## Usage

``` r
getParams(sim, time_range, geometric_mean = FALSE)

initialParams(sim)

finalParams(sim)
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
`initial_n_other`, and `initial_effort` set to the values from the
selected time of the simulation.

## Details

`getParams()` returns the state averaged over a chosen `time_range`, or
at a single time point. When no `time_range` is given, the state at the
final time step is returned.

`initialParams()` returns the state at the initial time of the
simulation, i.e., the `MizerParams` object that the simulation started
from.

`finalParams()` returns the state at the *last* saved time step. It is a
convenience wrapper around `getParams()` with no `time_range` argument.

The abundances set by `getParams()` are averages over the selected time
range. By default this is an arithmetic mean; set
`geometric_mean = TRUE` to use a geometric mean instead (this does not
affect the effort or other components, which are always averaged
arithmetically).

## Examples

``` r
sim <- project(NS_params, t_max = 20, effort = 0.5)
#> ℹ No `a` column so using a = 0.01 in w = a l^b, with w in g and l in cm.
#> ℹ No `b` column so using the isometric default b = 3 in w = a l^b.
# State at a specific time
params_10 <- getParams(sim, time_range = 10)
# State averaged over the last 10 years
params_avg <- getParams(sim, time_range = c(10, 20))
# State at the start and at the end of the simulation
params_start <- initialParams(sim)
params_end <- finalParams(sim)
```
