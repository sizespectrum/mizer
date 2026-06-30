# Extract the projection parameters used to produce a simulation

Returns the named list of arguments passed to
[`project()`](https://sizespectrum.org/mizer/reference/project.md) or
[`projectToSteady()`](https://sizespectrum.org/mizer/reference/projectToSteady.md)
when producing this `MizerSim` object, such as `method` and `dt`.
Returns an empty list for simulations produced by older versions of
mizer.

## Usage

``` r
getSimParams(sim)
```

## Arguments

- sim:

  A MizerSim object

## Value

A named list of projection parameters.

## Examples

``` r
sim <- project(NS_params, t_max = 0.1, dt = 0.05, method = "predictor-corrector")
getSimParams(sim)
#> $method
#> [1] "predictor_corrector"
#> 
#> $dt
#> [1] 0.05
#> 
#> $flux_limiter
#> [1] "none"
#> 
```
