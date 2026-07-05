# Calculate a selected subset of the rates

Internal helper used by the `MizerSim` rate getters. Given rate
functions already resolved once with `projectRateFunctions()`, it
calculates only those rates needed to obtain the requested `targets`
(plus their dependencies), avoiding both the per-time-step cost of
re-resolving the functions and the cost of computing rates that are not
required. The individual calculations mirror those in
[`mizerRates()`](https://sizespectrum.org/mizer/reference/mizerRates.md).

## Usage

``` r
mizer_rates_subset(
  params,
  n,
  n_pp,
  n_other,
  t,
  effort,
  rates_fns,
  targets,
  ...
)
```

## Arguments

- params:

  A valid `MizerParams` object.

- n:

  A matrix of species abundances (species x size).

- n_pp:

  A vector of the resource abundance by size.

- n_other:

  A named list of the abundances of other components.

- t:

  The time for the calculation.

- effort:

  The fishing effort. Only used when a target requires the fishing
  mortality.

- rates_fns:

  Named list of resolved rate functions, as returned by
  `projectRateFunctions()`.

- targets:

  Character vector of rate names (as in `params@rates_funcs`) to
  calculate.

- ...:

  Passed on to the individual rate functions.

## Value

A named list of the calculated rates, using the same element names as
the list returned by
[`mizerRates()`](https://sizespectrum.org/mizer/reference/mizerRates.md).
