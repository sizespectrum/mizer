# Collect the extent of each labelled dimension

Internal helper for
[`sizeIntegral()`](https://sizespectrum.org/mizer/reference/sizeIntegral.md)
that checks that arrays sharing a dimension label agree on its extent.

## Usage

``` r
dim_extents(arrays, labels)
```

## Arguments

- arrays:

  A list of arrays.

- labels:

  A list of the corresponding label vectors.

## Value

A named integer vector giving the extent of each label.
