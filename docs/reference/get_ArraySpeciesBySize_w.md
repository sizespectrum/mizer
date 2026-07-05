# Get the size grid for an ArraySpeciesBySize object

Internal helper that returns the consumer size grid `params@w` or the
full prey/resource size grid `params@w_full`, depending on the number of
columns in the array.

## Usage

``` r
get_ArraySpeciesBySize_w(x)
```

## Arguments

- x:

  An `ArraySpeciesBySize` object.

## Value

A numeric vector giving the size represented by each column. When the
array is tagged as a bin average (`representation = "average"`) *and*
the model uses second-order bin-averaging
(`second_order_w[["bin_average"]]`), the geometric bin centres are
returned instead of the left bin edges, so that bin-averaged quantities
are drawn at the size where they actually live (see
[`bin_midpoints()`](https://sizespectrum.org/mizer/reference/bin_midpoints.md)).
Point-valued quantities and first-order models are unaffected, keeping
default plots unchanged.
