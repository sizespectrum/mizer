# Make a plot of the relative difference between two data frames

Used internally by
[`plotSpectraRelative()`](https://sizespectrum.org/mizer/reference/plotSpectraRelative.md)
and similar functions. The two data frames are joined on their shared
variables and the relative difference of their y-values is plotted
against the x-variable.

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
  size_axis = NULL
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

- size_axis:

  Optional. If non-NULL, the x-axis is converted to weight (`"w"`) or
  length (`"l"`).

## Value

A `mizer_plot` (ggplot2) object showing the relative difference.
