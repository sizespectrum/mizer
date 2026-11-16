# Check whether a model is at steady state

**\[experimental\]** Returns the verdict from
[`getSteadyAssessment()`](https://sizespectrum.org/mizer/reference/getSteadyAssessment.md):
`TRUE` if all required owners pass, `FALSE` if any measure fails, or
`NA` if no measure fails but required convergence support is missing.
The default assesses the whole model, including extension components.
Callback errors propagate.

## Usage

``` r
isSteady(
  params,
  tol = 0.05,
  effort = params@initial_effort,
  scope = c("whole", "conditional"),
  fixed = character(),
  criterion = c("state", "biomass"),
  tolerances = NULL,
  atol = 1e-12,
  t = 0,
  dt = 1e-04,
  ...
)
```

## Arguments

- params:

  A MizerParams object or an extension thereof.

- tol:

  Relative tolerance for built-in checks, in 1/year.

- effort:

  Fishing effort at which to evaluate steadiness.

- scope:

  `"whole"` or `"conditional"`.

- fixed:

  Canonical owner IDs held fixed under conditional analysis:
  `"consumers"`, `"resource"`, or `"component:name"`.

- criterion:

  `"state"` for absolute state residuals, or `"biomass"` for net biomass
  drift in built-in checks.

- tolerances:

  Optional named vector of owner/measure tolerance overrides.

- atol:

  Absolute tolerance for built-in checks, in grams per year.

- t:

  Time at which to evaluate the model, in years.

- dt:

  Step length for resource dynamics supplied as a step map.

- ...:

  Additional arguments made available to convergence callbacks.

## Value

A single logical value, possibly `NA`. Conditional verdicts carry
`scope` and `fixed` attributes.

## Details

The default state-residual criterion prevents changes at different sizes
from cancelling. Net biomass drift remains available explicitly with
`criterion = "biomass"`. A conditional assessment requires an explicit
list of owners held fixed; its scalar verdict retains `scope` and
`fixed` attributes so that it cannot be mistaken for a whole-model
judgement.

## See also

[`getSteadyAssessment()`](https://sizespectrum.org/mizer/reference/getSteadyAssessment.md),
[`getSteadyResidual()`](https://sizespectrum.org/mizer/reference/getSteadyResidual.md),
[`tuneSteadyState()`](https://sizespectrum.org/mizer/reference/tuneSteadyState.md),
[`findSteadyState()`](https://sizespectrum.org/mizer/reference/findSteadyState.md)

## Examples

``` r
isSteady(NS_params)
#> [1] TRUE
```
