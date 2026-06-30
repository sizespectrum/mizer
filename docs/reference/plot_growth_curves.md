# Build the growth-curves plot

Internal worker shared by the
[`plotGrowthCurves()`](https://sizespectrum.org/mizer/reference/plotGrowthCurves.md)
methods. It computes the modelled size at age, optionally adds a von
Bertalanffy curve and observed size-at-age data, and draws the plot.

## Usage

``` r
plot_growth_curves(
  params,
  species,
  max_age,
  percentage,
  species_panel,
  highlight,
  log_x,
  log_y,
  size_at_age,
  return_data
)
```

## Arguments

- params:

  A MizerParams object.

- species:

  The species to be plotted.

- max_age:

  The age up to which to plot the growth curve.

- percentage:

  If `TRUE`, size is shown as a percentage of maximum size.

- species_panel:

  If `TRUE`, each species is shown in its own panel.

- highlight:

  Name or vector of names of species to be highlighted.

- log_x, log_y:

  Logical flags for log10 axes.

- size_at_age:

  Optional data frame of observed size-at-age data.

- return_data:

  If `TRUE`, return the plotting data frame instead of the plot.

## Value

A `mizer_plot` (ggplot2) object, or the plotting data frame if
`return_data = TRUE`.
