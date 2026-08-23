# The complete plotting data of a resource-by-size array

The resource analogue of
[`ArraySpeciesBySize_plot_data()`](https://sizespectrum.org/mizer/reference/ArraySpeciesBySize_plot_data.md):
the weight limits, the conversion of the values and of the size
coordinate onto the requested axis, and the length limits, all done with
the array's own `params`. A resource array holds a single spectrum, so
there is no selection, no background and no total to form.

## Usage

``` r
ArrayResourceBySize_plot_data(
  x,
  wlim = c(NA, NA),
  llim = c(NA, NA),
  size_axis = "w",
  per_log_size = NULL
)
```

## Arguments

- x:

  An `ArrayResourceBySize` object.

- wlim:

  Numeric vector of length two giving the weight limits.

- llim:

  Numeric vector of length two giving the length limits, applied only on
  a length axis.

- size_axis:

  Either `"w"` (weight) or `"l"` (length).

- per_log_size:

  Whether to express a density per logarithmic size.

## Value

A data frame with the size coordinate in its first column, the values in
its second, and `Species` and `Legend` columns.
