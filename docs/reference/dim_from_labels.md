# The dimensions of an array, given its labels

Internal helper for
[`sizeIntegral()`](https://sizespectrum.org/mizer/reference/sizeIntegral.md).
A scalar has no labels and no dimensions, a vector has one.

## Usage

``` r
dim_from_labels(x, labels)
```

## Arguments

- x:

  An array, vector or scalar.

- labels:

  Its dimension labels.

## Value

An integer vector of dimensions, of the same length as `labels`.
