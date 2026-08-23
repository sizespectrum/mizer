# Make a plot comparing two data frames

Used internally by the comparison plotting functions such as
[`plotSpectra2()`](https://sizespectrum.org/mizer/reference/plotSpectra2.md)
and
[`plotCDF2()`](https://sizespectrum.org/mizer/reference/plotCDF2.md).
The two data frames are combined and drawn with colour identifying the
species or group and linetype identifying the object.

## Usage

``` r
plotComparisonDataFrame(
  frame1,
  frame2,
  params,
  name1 = "First",
  name2 = "Second",
  xlab = waiver(),
  ylab = waiver(),
  xtrans = "identity",
  ytrans = "identity",
  xlim = c(NA, NA),
  ylim = c(NA, NA),
  y_ticks = 6,
  highlight = NULL,
  legend_var = "Legend"
)
```

## Arguments

- frame1, frame2:

  Data frames sharing the same first three variables (x, y and grouping
  variable). The names of `frame1` are imposed on `frame2`.

- params:

  A MizerParams object, used for the line colours.

- name1, name2:

  Labels for the two data frames, used in the linetype legend.

- xlab, ylab:

  Labels for the x and y axes.

- xtrans, ytrans:

  Transformations for the x and y axes, e.g. `"log10"` or `"identity"`.

- xlim, ylim:

  Numeric vectors of length two giving the axis limits. Use `NA` to
  refer to the existing minimum or maximum.

- y_ticks:

  The approximate number of ticks desired on the y axis.

- highlight:

  Name or vector of names of the species to be highlighted.

- legend_var:

  Name of the variable used in the legend and to determine the line
  colour.

## Value

A `mizer_plot` (ggplot2) object.

## Details

Both data frames must arrive ready to plot: on the axis they will be
drawn against, with any total line already among their rows. Each
operand is prepared by whatever produced it, using *its own* model,
because a length axis and a density Jacobian both depend on the
weight-length relationship of the model the values came from. Doing it
here instead would silently impose the first model's parameters on the
second.
