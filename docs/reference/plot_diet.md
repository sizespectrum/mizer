# Build the diet-composition plot

Internal worker shared by the
[`plotDiet()`](https://sizespectrum.org/mizer/reference/plotDiet.md)
methods. It melts the diet array into a data frame, restricts to
meaningful size ranges, converts to a length axis if requested and draws
a stacked-area plot of prey proportions.

## Usage

``` r
plot_diet(
  params,
  n,
  diet,
  species,
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

- n:

  Array of species abundances (species by size).

- diet:

  Array of diet proportions (predator by size by prey).

- species:

  The predator species to be plotted.

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
