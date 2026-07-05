# Get the size grid for an ArrayTimeBySpeciesBySize object

Internal helper, the three-dimensional analogue of
[`get_ArraySpeciesBySize_w()`](https://sizespectrum.org/mizer/reference/get_ArraySpeciesBySize_w.md).
Returns the geometric bin centres (see
[`bin_midpoints()`](https://sizespectrum.org/mizer/reference/bin_midpoints.md))
when the array is tagged as a bin average and the model uses
second-order bin-averaging, otherwise the grid nodes read from the size
dimension names. Falls back to the dimension names when no `params` is
attached.

## Usage

``` r
get_ArrayTimeBySpeciesBySize_w(x)
```

## Arguments

- x:

  An `ArrayTimeBySpeciesBySize` object.

## Value

A numeric vector giving the size represented by each size slice.
