# Convert the contributors to a total and append the total line

Wraps the two steps that every plot of an array with `total = TRUE`
needs: the contributors are put on the same size axis as the data they
will join, and then summed there, see
[`add_total_line()`](https://sizespectrum.org/mizer/reference/add_total_line.md).

## Usage

``` r
append_total_line(
  plot_dat,
  total_dat,
  params,
  size_axis,
  x,
  per_log_size = NULL
)
```

## Arguments

- plot_dat:

  The converted plotting data to append the total to.

- total_dat:

  The unconverted contributors, see
  [`total_contributors()`](https://sizespectrum.org/mizer/reference/total_contributors.md).

- params:

  A MizerParams object.

- size_axis:

  Either `"w"` (weight) or `"l"` (length).

- x:

  The mizer array the data came from, which says whether the values are
  a density.

- per_log_size:

  Whether to express a density per logarithmic size.

## Value

`plot_dat` with the total appended as a series named `"Total"`.
