# Determine which rates must be calculated to obtain a set of target rates

Internal helper returning the transitive closure of `targets` over
`.rate_dependencies`, in an order in which the rates can be calculated
(each rate appears after all the rates it depends on).

## Usage

``` r
needed_rates(targets)
```

## Arguments

- targets:

  Character vector of rate names (as in `params@rates_funcs`).

## Value

A character vector of rate names.
