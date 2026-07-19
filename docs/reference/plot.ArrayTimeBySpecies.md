# Plot method for `ArrayTimeBySpecies` objects

See [`plot()`](https://sizespectrum.org/mizer/reference/plot.md) for an
overview of the mizer plotting system and the arguments shared by all of
its methods.

## Usage

``` r
# S3 method for class 'ArrayTimeBySpecies'
plot(
  x,
  species = NULL,
  tlim = c(NA, NA),
  y_ticks = 6,
  ylim = c(NA, NA),
  total = FALSE,
  background = TRUE,
  highlight = NULL,
  log_x = FALSE,
  log_y = TRUE,
  log = NULL,
  return_data = FALSE,
  ...
)
```

## Arguments

- x:

  An `ArrayTimeBySpecies` object.

- species:

  Character vector of species to include. `NULL` (default) means all
  species.

- tlim:

  A numeric vector of length two providing lower and upper limits for
  the time axis, e.g. `c(1980, 2000)`. Use `NA` to apply no limit at
  that end. Default is `c(NA, NA)`.

- y_ticks:

  The approximate number of ticks desired on the y axis.

- ylim:

  A numeric vector of length two providing lower and upper limits for
  the value (y) axis. Use `NA` to refer to the existing minimum or
  maximum.

- total:

  A boolean value that determines whether the total over all selected
  species is plotted as well. Default is `FALSE`.

- background:

  A boolean value that determines whether background species are
  included. Ignored if the model does not contain background species.
  Default is `TRUE`.

- highlight:

  Name or vector of names of the species to be highlighted.

- log_x:

  If `TRUE`, use a log10 x-axis. Default is `FALSE`.

- log_y:

  If `TRUE`, use a log10 y-axis. Default is `TRUE`.

- log:

  Character string specifying which axes should use log10 scales, in the
  same form as the base
  [`plot()`](https://sizespectrum.org/mizer/reference/plot.md) argument.
  For example, `"x"`, `"y"`, `"xy"` or `""`. If supplied, this overrides
  `log_x` and `log_y`.

- return_data:

  If `TRUE`, return the data frame instead of the plot.

- ...:

  Unused.

## Value

A ggplot2 object, unless `return_data = TRUE`, in which case a data
frame is returned.

## Examples

``` r
# \donttest{
plot(getBiomass(NS_sim))

plot(getBiomass(NS_sim), species = c("Cod", "Herring"), total = TRUE)

plot(getYield(NS_sim), species = c("Cod", "Herring"))

# }
```
