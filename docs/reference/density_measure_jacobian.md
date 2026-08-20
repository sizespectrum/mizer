# Jacobian converting between two density measures

Jacobian converting between two density measures

## Usage

``` r
density_measure_jacobian(from, to, w, l, b)
```

## Arguments

- from, to:

  Density measures, see
  [density_measures](https://sizespectrum.org/mizer/reference/density_measures.md).

- w, l, b:

  Numeric vectors of the same length giving the weight, the
  corresponding length, and the exponent of the weight-length
  relationship.

## Value

A numeric vector by which to multiply a density with respect to `from`
to obtain the density with respect to `to`.
