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
  llim = c(NA, NA),
  ylim = c(NA, NA),
  size_axis = c("w", "l"),
  per_log_size = NULL,
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
  using the weight-length relationship in
  [`resource_params()`](https://sizespectrum.org/mizer/reference/resource_params.md).

- per_log_size:

  For an array that holds a density, whether to plot it per logarithmic
  size (`TRUE`) rather than per size (`FALSE`). The default, `NULL`,
  plots the density as it stands. An error for an array that does not
  hold a density.

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

# }
```
