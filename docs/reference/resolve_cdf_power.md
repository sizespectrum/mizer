# Resolve the power of weight for a cumulative distribution

As
[`resolve_spectrum_power()`](https://sizespectrum.org/mizer/reference/resolve_spectrum_power.md),
except that `per_log_size = TRUE` is rejected: integrating a density
with respect to logarithmic size gives the same cumulative quantity as
integrating the corresponding density with respect to size, so the flag
would be meaningless here.

## Usage

``` r
resolve_cdf_power(power = NULL, biomass = NULL, per_log_size = NULL)
```

## Arguments

- power:

  The power of weight multiplying the number density, or `NULL`.

- biomass:

  Whether to plot a biomass density rather than a number density, or
  `NULL`. The default is `TRUE`.

- per_log_size:

  Whether to plot a density with respect to logarithmic size rather than
  with respect to size, or `NULL`. The default is `FALSE`.

## Value

A list with entries `power`, `biomass` and `per_log_size`.
