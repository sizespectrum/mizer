# Check whether a model is at steady state

**\[experimental\]** Returns `TRUE` if the model is at its steady state
(within a specified tolerance), `FALSE` otherwise.

## Usage

``` r
isSteady(params, tol = 0.05, effort = params@initial_effort, ...)
```

## Arguments

- params:

  A
  [MizerParams](https://sizespectrum.org/mizer/reference/MizerParams-class.md)
  object or an extension thereof.

- tol:

  Tolerance for the relative rate of biomass change in 1/year. Defaults
  to `0.05` (5% change per year).

- effort:

  The fishing effort at which to evaluate steadiness. By default the
  initial effort stored in `params`.

- ...:

  Additional arguments passed to methods.

## Value

`TRUE` if the model's biomass drift is within `tol`, `FALSE` otherwise.

## Details

Steadiness is judged by computing the relative rate of change of biomass
across all consumer species, resource, and other components (see
[`getSteadyResidual()`](https://sizespectrum.org/mizer/reference/getSteadyResidual.md)).
If the largest biomass drift is less than or equal to `tol`, the model
is considered to be at steady state.

## See also

[`getSteadyResidual()`](https://sizespectrum.org/mizer/reference/getSteadyResidual.md),
[`steady()`](https://sizespectrum.org/mizer/reference/steady.md),
[`steadyNewton()`](https://sizespectrum.org/mizer/reference/steadyNewton.md)

## Examples

``` r
isSteady(NS_params)
#> [1] TRUE

# \donttest{
# Moving a species abundance off its steady state makes isSteady() FALSE
params <- NS_params
initialN(params)[1, ] <- initialN(params)[1, ] * 2
isSteady(params)
#> [1] FALSE
# }
```
