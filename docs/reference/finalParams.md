# Extract the final state from a simulation

Returns the `MizerParams` object underlying the simulation with its
initial abundances set to the abundances at the *last* saved time step
of the simulation. This is a convenience wrapper around
[`getParams()`](https://sizespectrum.org/mizer/reference/getParams.md)
with no `time_range` argument (the default).

## Usage

``` r
finalParams(sim)
```

## Arguments

- sim:

  A `MizerSim` object.

## Value

A `MizerParams` object with initial values taken from the final time
step of the simulation.

## See also

[`getParams()`](https://sizespectrum.org/mizer/reference/getParams.md),
[`initialParams()`](https://sizespectrum.org/mizer/reference/initialParams.md)

## Examples

``` r
sim <- project(NS_params, t_max = 20, effort = 0.5)
params_end <- finalParams(sim)
```
