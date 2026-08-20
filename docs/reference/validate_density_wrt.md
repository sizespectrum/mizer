# Validate a density measure

Validate a density measure

## Usage

``` r
validate_density_wrt(density_wrt)
```

## Arguments

- density_wrt:

  A density measure, see
  [density_measures](https://sizespectrum.org/mizer/reference/density_measures.md).
  `NULL` and `NA` both stand for "not a density".

## Value

The validated measure, or `NA_character_` when the values are not a
density.
