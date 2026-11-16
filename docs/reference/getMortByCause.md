# Diagnose mortality by cause

Calculate mortality contributions separately from the authoritative
total mortality used in the dynamical equations. The getter evaluates
the complete effective rate pipeline once at each requested state, then
checks a declared complete decomposition against its total mortality.

## Usage

``` r
getMortByCause(
  object,
  rtol = 1e-08,
  atol = 1e-12,
  ...,
  predation_by_species = FALSE
)
```

## Arguments

- object:

  A MizerParams or MizerSim object.

- rtol:

  Relative tolerance for the decomposition check. Default `1e-8`.

- atol:

  Absolute tolerance in 1/year. Default `1e-12`.

- ...:

  For MizerParams, `n`, `n_pp`, `n_other` and `effort` default to the
  model's initial values and `t` defaults to zero. For MizerSim,
  `time_range` selects saved times (all by default), `t` is an alias for
  `time_range`, and `drop = TRUE` drops dimensions of length one, as in
  [`getMort()`](https://sizespectrum.org/mizer/reference/getMort.md).
  Simulation states and effort are taken from each saved time. Other
  arguments are passed to rate and diagnostic callbacks.

- predation_by_species:

  Whether to refine built-in predation into one cause per predator
  species. The refinement is used only when its sum agrees with the
  effective predation mortality at every species and size; otherwise the
  aggregate predation cause is retained. For the built-in callback and
  its public alias, any selected simulation time requiring aggregate
  predation makes every selected time use that grouping. Custom
  callbacks are responsible for stable grouping across times. Default
  `FALSE`.

## Value

A `MizerMortalityBreakdown` list with `total`, a named list of `causes`,
plotting `metadata`, logical `complete`, an `unavailable` data frame
including time and reasons, and `residual` (total minus the sum of
supplied causes). `total`, `causes` and `residual` use the same shapes
as [`getMort()`](https://sizespectrum.org/mizer/reference/getMort.md).
The residual is a diagnostic difference, not an attributed cause.
Simulation causes and metadata must remain consistent through time.

## Details

Built-in mortality is decomposed into predation, fishing, external
mortality and each additive contribution registered with
[`setComponent()`](https://sizespectrum.org/mizer/reference/setComponent.md)
or
[`other_mort()`](https://sizespectrum.org/mizer/reference/other_mort.md).
Replacing the outer `Rates` or `Mort` provider requires an explicitly
adopted decomposition through
`setProcess(params, "consumers", mortality_fun = "myBreakdown")`.
Unsupported decomposition returns an incomplete report with a reason and
the authoritative total; it does not assign unexplained mortality to a
cause.

A `mortality_fun` receives the model, `n`, `n_pp`, `n_other`, `effort`,
`t`, the complete `rates` list, `total` and `process = "consumers"`. It
returns a list with `causes` (a named list of finite numeric
species-by-size matrices) and `complete` (one logical value). Optional
`metadata` is a data frame with `id` and any of `label`, `owner`,
`colour`, `linetype`; `id` must name each cause once. Optional
`unavailable` is a data frame with `id`, `owner` and `reason` columns.
It is only allowed for incomplete decompositions. Callback errors and
malformed results are errors. A complete decomposition must satisfy
`abs(total - sum(causes)) <= atol + rtol * abs(total)` at every species,
size and saved time.

## See also

[`getMort()`](https://sizespectrum.org/mizer/reference/getMort.md),
[`setProcess()`](https://sizespectrum.org/mizer/reference/setProcess.md)

Other summary functions:
[`getBiomass()`](https://sizespectrum.org/mizer/reference/getBiomass.md),
[`getDiet()`](https://sizespectrum.org/mizer/reference/getDiet.md),
[`getEncounterByPrey()`](https://sizespectrum.org/mizer/reference/getEncounterByPrey.md),
[`getGrowthCurves()`](https://sizespectrum.org/mizer/reference/getGrowthCurves.md),
[`getN()`](https://sizespectrum.org/mizer/reference/getN.md),
[`getSSB()`](https://sizespectrum.org/mizer/reference/getSSB.md),
[`getSteadyResidual()`](https://sizespectrum.org/mizer/reference/getSteadyResidual.md),
[`getTrophicLevel()`](https://sizespectrum.org/mizer/reference/getTrophicLevel.md),
[`getTrophicLevelBySpecies()`](https://sizespectrum.org/mizer/reference/getTrophicLevelBySpecies.md),
[`getYield()`](https://sizespectrum.org/mizer/reference/getYield.md),
[`getYieldGear()`](https://sizespectrum.org/mizer/reference/getYieldGear.md)

## Examples

``` r
mortality <- getMortByCause(NS_params)
mortality$complete
#> [1] TRUE
mortality$metadata
#>          id     label     owner  colour linetype
#> 1 predation Predation consumers #D55E00    solid
#> 2   fishing   Fishing consumers     red    solid
#> 3  external  External consumers    grey    solid
plot(mortality$total)
```
