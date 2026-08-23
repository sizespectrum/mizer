# The complete plotting data of a species-by-size array

Everything a plot of an `ArraySpeciesBySize` needs, prepared once: the
species selection, the background grouping, the masking of sizes outside
a species' own range, the weight limits, the conversion of the values
and of the size coordinate onto the requested axis, the total line, and
the length limits.

## Usage

``` r
ArraySpeciesBySize_plot_data(
  x,
  species = NULL,
  all.sizes = FALSE,
  wlim = c(NA, NA),
  llim = c(NA, NA),
  total = FALSE,
  background = TRUE,
  size_axis = "w",
  per_log_size = NULL
)
```

## Arguments

- x:

  An `ArraySpeciesBySize` object.

- species:

  Character vector of species to include, or `NULL` for all.

- all.sizes:

  If `FALSE`, values outside a species' size range are removed.

- wlim:

  Numeric vector of length two giving the weight limits.

- llim:

  Numeric vector of length two giving the length limits, applied only on
  a length axis.

- total:

  Whether to append the total line, see
  [`total_contributors()`](https://sizespectrum.org/mizer/reference/total_contributors.md).

- background:

  Whether background species are included.

- size_axis:

  Either `"w"` (weight) or `"l"` (length).

- per_log_size:

  Whether to express a density per logarithmic size.

## Value

A data frame with the size coordinate in its first column, the values in
its second, and `Species` and `Legend` columns.

## Details

All of it uses the array's *own* `params`. That matters for the
comparison plots, where the two operands may come from different models:
a length axis and a density Jacobian are both built from the
weight-length relationship of the model the values came from, so
preparing the second array with the first one's parameters would put it
in the wrong place on the axis. Each operand is therefore prepared here,
on its own, and the comparison renderers receive data that is already on
the axis it will be drawn against.
