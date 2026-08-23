# Make a plot of the relative difference between two data frames

Used internally by
[`plotSpectraRelative()`](https://sizespectrum.org/mizer/reference/plotSpectraRelative.md)
and similar functions. The two data frames are matched up on their
shared variables and the relative difference of their y-values is
plotted against the x-variable.

## Usage

``` r
plotRelativeDataFrame(
  frame1,
  frame2,
  params,
  xlab = waiver(),
  xtrans = "identity",
  xlim = c(NA, NA),
  ylim = c(NA, NA),
  highlight = NULL,
  legend_var = "Legend",
  interpolate = FALSE
)
```

## Arguments

- frame1, frame2:

  Data frames sharing the same first three variables (x, y and grouping
  variable). The names of `frame1` are imposed on `frame2`.

- params:

  A MizerParams object, used for the line colours.

- xlab:

  Label for the x-axis.

- xtrans:

  Transformation for the x-axis, e.g. `"log10"` or `"identity"`.

- xlim, ylim:

  Numeric vectors of length two giving the axis limits. Use `NA` to
  refer to the existing minimum or maximum.

- highlight:

  Name or vector of names of the species to be highlighted.

- legend_var:

  Name of the variable used in the legend and to determine the line
  colour.

- interpolate:

  Whether the two series may sit on different x-grids and should be
  interpolated onto a common one, see
  [`interpolate_relative_frames()`](https://sizespectrum.org/mizer/reference/interpolate_relative_frames.md).
  `TRUE` for a size axis, where the two models can convert weight to
  length differently; `FALSE` for a time axis, where they share the
  saved times or share nothing.

## Value

A `mizer_plot` (ggplot2) object showing the relative difference.

## Details

Both data frames must arrive ready to plot, on the axis they will be
drawn against and with any total line already among their rows. See
[`plotComparisonDataFrame()`](https://sizespectrum.org/mizer/reference/plotComparisonDataFrame.md)
for why the preparation belongs to whatever produced each operand rather
than here.
