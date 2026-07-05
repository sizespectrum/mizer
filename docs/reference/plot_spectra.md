# Build the size-spectrum plot

Internal worker shared by the
[`plotSpectra()`](https://sizespectrum.org/mizer/reference/plotSpectra.md)
methods. It assembles the plotting data frame from the species and
resource abundances, applies the size and abundance limits, optionally
converts to a length axis, and either returns the data or draws the plot
via
[`plotDataFrame()`](https://sizespectrum.org/mizer/reference/plotDataFrame.md).

## Usage

``` r
plot_spectra(
  params,
  n,
  n_pp,
  species,
  wlim,
  llim,
  ylim,
  power,
  total,
  resource,
  background,
  highlight,
  log_x,
  log_y,
  size_axis,
  return_data
)
```

## Arguments

- params:

  A MizerParams object.

- n:

  Array of species abundances (species by size).

- n_pp:

  Vector of resource abundance.

- species:

  The species to be plotted.

- wlim, llim:

  Numeric vectors of length two giving the weight and length limits.

- ylim:

  Numeric vector of length two giving the y-axis limits.

- power:

  The abundance is multiplied by weight raised to this power.

- total:

  Whether to include the total community abundance.

- resource:

  Whether to include the resource spectrum.

- background:

  Whether to include background species.

- highlight:

  Name or vector of names of species to be highlighted.

- log_x, log_y:

  Logical flags for log10 axes.

- size_axis:

  Either `"w"` (weight) or `"l"` (length).

- return_data:

  If `TRUE`, return the plotting data frame instead of the plot.

## Value

A `mizer_plot` (ggplot2) object, or the plotting data frame if
`return_data = TRUE`.
