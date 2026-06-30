# Build the cumulative-distribution plot

Internal worker shared by the
[`plotCDF()`](https://sizespectrum.org/mizer/reference/plotCDF.md)
methods. It integrates the spectra data over size (optionally
normalising), converts to a length axis if requested, and either returns
the data or draws the plot via
[`plotDataFrame()`](https://sizespectrum.org/mizer/reference/plotDataFrame.md).

## Usage

``` r
plot_cdf(
  plot_dat,
  params,
  power,
  normalise,
  log_x,
  log_y,
  wlim,
  llim,
  ylim,
  highlight,
  size_axis,
  return_data
)
```

## Arguments

- plot_dat:

  Spectra plotting data as produced for
  [`plotSpectra()`](https://sizespectrum.org/mizer/reference/plotSpectra.md).

- params:

  A MizerParams object.

- power:

  The power of weight that the abundance was multiplied by, used for the
  y-axis label.

- normalise:

  If `TRUE`, each curve is divided by its final value.

- log_x, log_y:

  Logical flags for log10 axes.

- wlim, llim:

  Numeric vectors of length two giving the weight and length limits.

- ylim:

  Numeric vector of length two giving the y-axis limits.

- highlight:

  Name or vector of names of species to be highlighted.

- size_axis:

  Either `"w"` (weight) or `"l"` (length).

- return_data:

  If `TRUE`, return the cumulative-distribution data frame instead of
  the plot.

## Value

A `mizer_plot` (ggplot2) object, or the data frame if
`return_data = TRUE`.
