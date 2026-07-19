# Plot method for `ArrayResourceBySize` objects

See [`plot()`](https://sizespectrum.org/mizer/reference/plot.md) for an
overview of the mizer plotting system and the arguments shared by all of
its methods.

## Usage

``` r
# S3 method for class 'ArrayResourceBySize'
plot(
  x,
  return_data = FALSE,
  log_x = TRUE,
  log_y = TRUE,
  log = NULL,
  wlim = c(NA, NA),
  ylim = c(NA, NA),
  y_ticks = 6,
  ...
)
```

## Arguments

- x:

  An `ArrayResourceBySize` object.

- return_data:

  If `TRUE`, return the data frame instead of the plot.

- log_x:

  If `TRUE`, use a log10 x-axis. Default is `TRUE`.

- log_y:

  If `TRUE`, use a log10 y-axis. Default is `TRUE`.

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

- ylim:

  A numeric vector of length two providing lower and upper limits for
  the value (y) axis. Use `NA` to refer to the existing minimum or
  maximum.

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
plot(getResourceMort(NS_params))

plot(initialNResource(NS_params))
#> Warning: log-10 transformation introduced infinite values.

# }
```
