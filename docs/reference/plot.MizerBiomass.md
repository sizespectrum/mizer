# Plot a biomass measurement at one state

Explicit resource or component inclusion in
[`getBiomass()`](https://sizespectrum.org/mizer/reference/getBiomass.md)
returns a named numeric `MizerBiomass` vector carrying its own series
labels, colours and coverage. This method draws those measurements
without discovering additional series from the model. Unavailable
measurements remain in the returned data and are explained in the plot
caption. A total requires complete coverage.

## Usage

``` r
# S3 method for class 'MizerBiomass'
plot(
  x,
  species = NULL,
  total = FALSE,
  log_y = TRUE,
  ylim = c(NA, NA),
  return_data = FALSE,
  ...
)
```

## Arguments

- x:

  A `MizerBiomass` vector returned by
  [`getBiomass()`](https://sizespectrum.org/mizer/reference/getBiomass.md).

- species:

  Series identifiers to draw; NULL draws every measured series.

- total:

  Add the total of every measured series, before selection.

- log_y:

  Use a logarithmic biomass axis.

- ylim:

  Numeric lower and upper biomass-axis limits; NA chooses them.

- return_data:

  Return the plotting data instead of a plot.

- ...:

  Currently unused.

## Value

A mizer plot, or a data frame when `return_data = TRUE`.
