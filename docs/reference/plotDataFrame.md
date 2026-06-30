# Make a plot from a data frame

This is used internally by most plotting functions.

## Usage

``` r
plotDataFrame(
  frame,
  params,
  style = "line",
  xlab = waiver(),
  ylab = waiver(),
  xtrans = "identity",
  ytrans = "identity",
  xlim = c(NA, NA),
  ylim = c(NA, NA),
  y_ticks = 6,
  highlight = NULL,
  legend_var = NULL,
  wrap_var = NULL,
  wrap_scale = NULL
)
```

## Arguments

- frame:

  A data frame with at least three variables. The first three variables
  are used, in that order, as:

  1.  Variable to be plotted on x-axis

  2.  Variable to be plotted on y-axis

  3.  Grouping variable

- params:

  A MizerParams object, which is used for the line colours and line
  types.

- style:

  The style of the plot. Available options are `"line"` for
  [`geom_line()`](https://ggplot2.tidyverse.org/reference/geom_path.html)
  and `"area"` for
  [`geom_area()`](https://ggplot2.tidyverse.org/reference/geom_ribbon.html).
  Default is `"line"`.

- xlab:

  Label for the x-axis

- ylab:

  Label for the y-axis

- xtrans:

  Transformation for the x-axis. Often "log10" may be useful instead of
  the default of "identity".

- ytrans:

  Transformation for the y-axis.

- xlim:

  A numeric vector of length two providing lower and upper limits for
  the x axis. Use NA to refer to the existing minimum or maximum.

- ylim:

  A numeric vector of length two providing lower and upper limits for
  the y axis. Use NA to refer to the existing minimum or maximum.

- y_ticks:

  The approximate number of ticks desired on the y axis

- highlight:

  Name or vector of names of the species to be highlighted.

- legend_var:

  The name of the variable that should be used in the legend and to
  determine the line style. If NULL then the grouping variable is used
  for this purpose.

- wrap_var:

  Optional. The name of the variable that should be used for creating
  wrapped facets.

- wrap_scale:

  Optional. Used to pass the scales argument to facet_wrap().

## Value

A ggplot2 object
