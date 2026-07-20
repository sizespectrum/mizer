# Plot method for `ArrayTimeByResourceBySize` objects

See [`plot()`](https://sizespectrum.org/mizer/reference/plot.md) for an
overview of the mizer plotting system. This method plots a single time
slice, by first extracting it as an `ArrayResourceBySize` object and
delegating to
[`plot.ArrayResourceBySize()`](https://sizespectrum.org/mizer/reference/plot.ArrayResourceBySize.md),
which the further arguments in `...` are passed on to.

## Usage

``` r
# S3 method for class 'ArrayTimeByResourceBySize'
plot(x, time = NULL, ...)
```

## Arguments

- x:

  An `ArrayTimeByResourceBySize` object.

- time:

  The time to display. Default (`NULL`) is the final time step.

- ...:

  Passed on to
  [`plot.ArrayResourceBySize()`](https://sizespectrum.org/mizer/reference/plot.ArrayResourceBySize.md).

## Value

A ggplot2 object, unless `return_data = TRUE`, in which case a data
frame is returned.

## Examples

``` r
# \donttest{
plot(NResource(NS_sim))
#> Warning: log-10 transformation introduced infinite values.

# }
```
