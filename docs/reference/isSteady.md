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
across all consumer species and the resource (see
[`getSteadyResidual()`](https://sizespectrum.org/mizer/reference/getSteadyResidual.md)).
If the largest biomass drift is less than or equal to `tol`, the model
is considered to be at steady state.

## Components are not included

Components registered with
[`setComponent()`](https://sizespectrum.org/mizer/reference/setComponent.md)
are deliberately left out of this judgement. Their state can be any
object at all, so mizer does not know what currency its entries are in
and cannot form a biomass for them; and
[`tuneSteadyState()`](https://sizespectrum.org/mizer/reference/tuneSteadyState.md),
[`findSteadyState()`](https://sizespectrum.org/mizer/reference/findSteadyState.md)
and
[`getStability()`](https://sizespectrum.org/mizer/reference/getStability.md)
all hold them fixed, so a criterion that included them could be one that
no mizer function is able to satisfy.

A model can therefore be `isSteady()` while a component of it is still
changing. If your model has components with dynamics of their own, check
them as well: `attr(getSteadyResidual(params), "other")` holds their
rates of change, and mizer names any that are moving whenever it reports
on the model's steady state. To settle them along with everything else,
project the model with
[`projectUntilSettled()`](https://sizespectrum.org/mizer/reference/projectUntilSettled.md),
which advances the components like every other state variable.

## See also

[`getSteadyResidual()`](https://sizespectrum.org/mizer/reference/getSteadyResidual.md),
[`tuneSteadyState()`](https://sizespectrum.org/mizer/reference/tuneSteadyState.md),
[`findSteadyState()`](https://sizespectrum.org/mizer/reference/findSteadyState.md)

## Examples

``` r
isSteady(NS_params)
#> ℹ No `a` column so using a = 0.01 in w = a l^b, with w in g and l in cm.
#> ℹ No `b` column so using the isometric default b = 3 in w = a l^b.
#> [1] TRUE

# \donttest{
# Moving a species abundance off its steady state makes isSteady() FALSE
params <- NS_params
initialN(params)[1, ] <- initialN(params)[1, ] * 2
#> ℹ No `a` column so using a = 0.01 in w = a l^b, with w in g and l in cm.
#> ℹ No `b` column so using the isometric default b = 3 in w = a l^b.
#> ℹ No `a` column so using a = 0.01 in w = a l^b, with w in g and l in cm.
#> ℹ No `b` column so using the isometric default b = 3 in w = a l^b.
isSteady(params)
#> ℹ No `a` column so using a = 0.01 in w = a l^b, with w in g and l in cm.
#> ℹ No `b` column so using the isometric default b = 3 in w = a l^b.
#> [1] FALSE
# }
```
