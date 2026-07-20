# Plot method for `ArraySpeciesBySize` objects

See [`plot()`](https://sizespectrum.org/mizer/reference/plot.md) for an
overview of the mizer plotting system and the arguments shared by all of
its methods.

## Usage

``` r
# S3 method for class 'ArraySpeciesBySize'
plot(
  x,
  species = NULL,
  all.sizes = FALSE,
  highlight = NULL,
  return_data = FALSE,
  log_x = TRUE,
  log_y = FALSE,
  log = NULL,
  wlim = c(NA, NA),
  llim = c(NA, NA),
  ylim = c(NA, NA),
  size_axis = c("w", "l"),
  total = FALSE,
  background = TRUE,
  y_ticks = 6,
  ...
)
```

## Arguments

- x:

  An `ArraySpeciesBySize` object.

- species:

  Character vector of species to include. `NULL` (default) means all
  species.

- all.sizes:

  If `FALSE` (default), values outside a species' size range (`w_min` to
  `w_max`) are removed.

- highlight:

  Name or vector of names of the species to be highlighted.

- return_data:

  If `TRUE`, return the data frame instead of the plot.

- log_x:

  If `TRUE`, use a log10 x-axis. Default is `TRUE`.

- log_y:

  If `TRUE`, use a log10 y-axis. Default is `FALSE`.

- log:

  Character string specifying which axes should use log10 scales, in the
  same form as the base
  [`plot()`](https://sizespectrum.org/mizer/reference/plot.md) argument.
  For example, `"x"`, `"y"`, `"xy"` or `""`. If supplied, this overrides
  `log_x` and `log_y`.

- wlim:

  A numeric vector of length two providing lower and upper limits for
  the weight (x) axis. Use `NA` to refer to the existing minimum or
  maximum.

- llim:

  A numeric vector of length two providing lower and upper limits for
  the length (x) axis when `size_axis = "l"`. Use `NA` to refer to the
  existing minimum or maximum.

- ylim:

  A numeric vector of length two providing lower and upper limits for
  the value (y) axis. Use `NA` to refer to the existing minimum or
  maximum.

- size_axis:

  Whether to plot size as weight (`"w"`, default) or length (`"l"`),
  using the allometric weight-length relationship.

- total:

  A boolean value that determines whether the total over all selected
  species is plotted as well. Default is `FALSE`.

- background:

  A boolean value that determines whether background species are
  included. Ignored if the model does not contain background species.
  Default is `TRUE`.

- y_ticks:

  The approximate number of ticks desired on the y axis.

- ...:

  Unused.

## Value

A ggplot2 object, unless `return_data = TRUE`, in which case a data
frame is returned.

## Examples

``` r
# \donttest{
plot(getEncounter(NS_params))

plot(getFeedingLevel(NS_params), species = c("Cod", "Herring"))

plot(getPredMort(NS_params), species = c("Cod", "Herring"),
     size_axis = "l")

# }
```
