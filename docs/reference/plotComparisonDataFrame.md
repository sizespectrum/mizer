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
  legend_var = "Legend",
  size_axis = NULL,
  density_wrt = NA_character_,
  per_log_size = NULL,
  total_dat = NULL
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

- size_axis:

  Optional. If non-NULL, the x-axis is converted to weight (`"w"`) or
  length (`"l"`).

- density_wrt:

  The measure the values are a density with respect to, see
  [density_measures](https://sizespectrum.org/mizer/reference/density_measures.md).
  `NA` (the default) means the values are not a density and are left
  alone when the size axis changes.

- per_log_size:

  Whether to express a density per logarithmic size. `NULL` (the
  default) keeps whichever the values already are.

- total_dat:

  Optional data frame of the contributors to a total, with a `Model`
  column identifying which of the two they belong to. The total is
  summed from them after the size axis has been converted, see
  [`add_total_line()`](https://sizespectrum.org/mizer/reference/add_total_line.md).

## Value

A `mizer_plot` (ggplot2) object.
