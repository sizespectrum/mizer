# Get the size grid for an ArrayTimeBySpeciesBySize object

Internal helper, the three-dimensional analogue of
[`get_ArraySpeciesBySize_w()`](https://sizespectrum.org/mizer/reference/get_ArraySpeciesBySize_w.md).
Returns the geometric bin centres (see
[`bin_midpoints()`](https://sizespectrum.org/mizer/reference/bin_midpoints.md))
when the array is tagged as a bin average and the model uses
second-order bin-averaging, otherwise the original grid nodes. The
stored size mapping preserves these coordinates through subsetting.
Falls back to dimension names when no model grid is available.

## Usage

``` r
get_ArrayTimeBySpeciesBySize_w(x)
```

## Arguments

- x:

  An `ArrayTimeBySpeciesBySize` object.

## Value

A numeric vector giving the size represented by each size slice.
