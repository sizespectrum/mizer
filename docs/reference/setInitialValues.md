# Set initial values to values from a simulation

**\[deprecated\]** This function is deprecated. Use
[`getParams()`](https://sizespectrum.org/mizer/reference/getParams.md),
[`initialParams()`](https://sizespectrum.org/mizer/reference/initialParams.md),
or
[`finalParams()`](https://sizespectrum.org/mizer/reference/finalParams.md)
instead. These functions return a `MizerParams` object with the
ecosystem state extracted from a simulation.

## Usage

``` r
setInitialValues(params, sim, time_range, geometric_mean = FALSE, ...)
```

## Arguments

- params:

  A `MizerParams` object in which to set the initial values

- sim:

  A `MizerSim` object from which to take the values.

- time_range:

  The time range to average the abundances over. Can be a vector of
  values, a vector of min and max time, or a single value. Only the
  range of times is relevant, i.e., all times between the smallest and
  largest will be selected. Default is the final time step.

- geometric_mean:

  **\[experimental\]** If TRUE then the average of the abundances over
  the time range is a geometric mean instead of the default arithmetic
  mean. This does not affect the average of the effort or of other
  components, which is always arithmetic.

- ...:

  Additional arguments passed to the method.

## Value

The `params` object with updated initial values and initial effort.

## Details

**\[deprecated\]**

## Examples

``` r
# \donttest{
params <- NS_params
sim <- project(params, t_max = 20, effort = 0.5)
params <- setInitialValues(params, sim)
#> Warning: `setInitialValues()` was deprecated in mizer 3.0.0.
#> ℹ Use `getParams(sim, time_range, geometric_mean)` to extract a MizerParams
#>   object with updated initial values. Convenience wrappers `initialParams()`
#>   and `finalParams()` extract the first and last time steps.
# }
```
