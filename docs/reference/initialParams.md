# Extract the initial state from a simulation

Returns the `MizerParams` object underlying the simulation with its
initial abundances set to the abundances at the initial time of the
simulation.

## Usage

``` r
initialParams(sim)
```

## Arguments

- sim:

  A `MizerSim` object.

## Value

A `MizerParams` object with initial state of the simulation.

## See also

[`getParams()`](https://sizespectrum.org/mizer/reference/getParams.md),
[`finalParams()`](https://sizespectrum.org/mizer/reference/finalParams.md)

## Examples

``` r
sim <- project(NS_params, t_max = 20, effort = 0.5)
params_start <- initialParams(sim)
```
