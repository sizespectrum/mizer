# Merge two ordered sets of dimension labels

Internal helper for
[`sizeIntegral()`](https://sizespectrum.org/mizer/reference/sizeIntegral.md).
Interleaves the labels of two arrays into the labels of the array
holding their product, keeping the relative order of the labels within
each of the two inputs. Labels that occur in both inputs occur once in
the result, which is how a weighting array with a `"time"` dimension is
lined up with the times of the abundance rather than multiplied out
against them.

## Usage

``` r
merge_dim_labels(a, b)
```

## Arguments

- a, b:

  Character vectors of dimension labels.

## Value

A character vector of labels containing each label of `a` and `b` once.
