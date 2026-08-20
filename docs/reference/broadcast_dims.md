# Expand an array to a larger set of labelled dimensions

Internal helper for
[`sizeIntegral()`](https://sizespectrum.org/mizer/reference/sizeIntegral.md).
Replicates `x` along the dimensions it does not have and permutes its
dimensions into the order given by `target_labels`, so that arrays with
different dimensions can be multiplied together elementwise.

## Usage

``` r
broadcast_dims(x, labels, target_labels, extent)
```

## Arguments

- x:

  An array, vector or scalar.

- labels:

  The dimension labels of `x`.

- target_labels:

  The dimension labels of the result. Must contain all of `labels`.

- extent:

  A named vector giving the extent of each label.

## Value

An array with dimensions `extent[target_labels]`.
