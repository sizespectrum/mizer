# Build the standard mortality diagnostic from supplied rates

An author helper for a registered `mortality_fun` callback. It
partitions the supplied effective rates into predation, fishing,
external mortality and registered additive contributions without
evaluating the rate pipeline again. An extension that changes the outer
mortality assembler can call this helper, append its additional terms,
and return the resulting specification.

## Usage

``` r
mizerMortalityBreakdown(
  params,
  n,
  n_pp,
  n_other,
  rates,
  t,
  effort,
  total,
  process,
  predation_by_species = FALSE,
  rtol = 1e-08,
  atol = 1e-12,
  ...
)
```

## Arguments

- params:

  A MizerParams object.

- n, n_pp, n_other:

  Current consumer, resource and component states.

- rates:

  The complete effective rate list supplied to the callback.

- t:

  Current model time.

- effort:

  Current fishing effort.

- total:

  Effective total mortality supplied to the callback.

- process:

  The diagnostic owner, `"consumers"`.

- predation_by_species:

  Whether to refine built-in predation into one cause per predator
  species. The refinement is used only when its sum agrees with the
  effective predation mortality at every species and size; otherwise the
  aggregate predation cause is retained. For the built-in callback and
  its public alias, any selected simulation time requiring aggregate
  predation makes every selected time use that grouping. Custom
  callbacks are responsible for stable grouping across times. Default
  `FALSE`.

- rtol:

  Relative tolerance for the decomposition check. Default `1e-8`.

- atol:

  Absolute tolerance in 1/year. Default `1e-12`.

- ...:

  Further callback context.

## Value

A raw callback specification with `causes`, `metadata` and
`complete = TRUE`. Authors must append every additional mortality term
before returning it as a complete decomposition.

## Details

This helper returns an unchecked callback specification. The complete
sum is checked against the effective total by
[`getMortByCause()`](https://sizespectrum.org/mizer/reference/getMortByCause.md)
after the author's callback returns, so terms added by that callback
participate in the check.

## See also

[`setProcess()`](https://sizespectrum.org/mizer/reference/setProcess.md),
[`getMortByCause()`](https://sizespectrum.org/mizer/reference/getMortByCause.md)
