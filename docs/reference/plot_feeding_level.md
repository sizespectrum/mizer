# Build the feeding-level plot

Internal worker shared by the
[`plotFeedingLevel()`](https://sizespectrum.org/mizer/reference/plotFeedingLevel.md)
methods. It assembles the plotting data, optionally adds the critical
feeding level, restricts to each species' size range, converts to a
length axis if requested and draws the plot.

## Usage

``` r
plot_feeding_level(
  params,
  feed,
  species,
  highlight,
  all.sizes,
  include_critical,
  log_x,
  log_y,
  wlim,
  llim,
  size_axis,
  return_data
)
```

## Arguments

- params:

  A MizerParams object.

- feed:

  Array of feeding levels (species by size).

- species:

  The species to be plotted.

- highlight:

  Name or vector of names of species to be highlighted.

- all.sizes:

  If `FALSE`, feeding levels outside each species' size range are
  removed.

- include_critical:

  Whether to also plot the critical feeding level.

- log_x, log_y:

  Logical flags for log10 axes.

- wlim, llim:

  Numeric vectors of length two giving the weight and length limits.

- size_axis:

  Either `"w"` (weight) or `"l"` (length).

- return_data:

  If `TRUE`, return the plotting data frame instead of the plot.

## Value

A `mizer_plot` (ggplot2) object, or the plotting data frame if
`return_data = TRUE`.
