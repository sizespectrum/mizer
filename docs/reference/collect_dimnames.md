# Assemble the dimnames of a broadcast array

Internal helper for
[`sizeIntegral()`](https://sizespectrum.org/mizer/reference/sizeIntegral.md)
that takes the dimnames of each labelled dimension from the first of the
given arrays that has them.

## Usage

``` r
collect_dimnames(target_labels, arrays, labels)
```

## Arguments

- target_labels:

  The dimension labels of the result.

- arrays:

  A list of arrays.

- labels:

  A list of the corresponding label vectors.

## Value

A named list of dimnames, or `NULL` if none of the arrays has any.
