# Plot method for `MizerScan` objects

**\[experimental\]** Draws the result of a
[`scanModel()`](https://sizespectrum.org/mizer/reference/scanModel.md)
run: the measured quantity against the quantity that was scanned, with a
band showing the range the quantity takes over the attractor wherever
that attractor is not a fixed point.

## Usage

``` r
# S3 method for class 'MizerScan'
plot(
  x,
  species = NULL,
  style = c("ribbon", "envelope", "line"),
  highlight = NULL,
  log_x = FALSE,
  log_y = TRUE,
  log = NULL,
  xlim = c(NA, NA),
  ylim = c(NA, NA),
  y_ticks = 6,
  reference_lines = TRUE,
  mark_max = FALSE,
  show_unsettled = TRUE,
  return_data = FALSE,
  ...
)
```

## Arguments

- x:

  A `MizerScan` object.

- species:

  The species to show. By default all series in the scan.

- style:

  One of `"ribbon"` (default: the average as a line inside the band),
  `"envelope"` (lines along the edges of the band, no average) or
  `"line"` (no band).

- highlight:

  Name or vector of names of the species to be highlighted with a
  thicker line.

- log_x, log_y, log:

  Whether to use logarithmic axes, see
  [`parsePlotLog()`](https://sizespectrum.org/mizer/reference/parsePlotLog.md).

- xlim, ylim:

  Numeric vectors of length two giving the axis limits. Use `NA` to
  refer to the existing minimum or maximum.

- y_ticks:

  The approximate number of ticks desired on the y axis.

- reference_lines:

  Whether to draw the reference lines stored in the scan, or a named
  numeric vector of x positions to draw instead.

- mark_max:

  Whether to mark, for each series, the scanned value at which the
  measured quantity is largest. See
  [`MizerScan()`](https://sizespectrum.org/mizer/reference/MizerScan.md).

- show_unsettled:

  Whether to mark the scan values where the model did not settle onto an
  attractor.

- return_data:

  Whether to return the data frame used for the plot instead of the plot
  itself.

- ...:

  Unused.

## Value

A ggplot2 object, unless `return_data = TRUE`, in which case the data
frame used for the plot is returned.

## Details

A model that settles on a fixed point contributes a single value, so the
band has zero width there. A model that settles on a limit cycle
contributes the average over one period as the line and the range over
that period as the band, so a Hopf bifurcation shows up as the scan
value at which the band opens up.

Scan values where the model reached neither a fixed point nor a limit
cycle within the time allowed are marked with a cross, because the value
plotted there is only an average over the last few years of a run that
was still changing.

## See also

[`scanModel()`](https://sizespectrum.org/mizer/reference/scanModel.md),
[`MizerScan()`](https://sizespectrum.org/mizer/reference/MizerScan.md),
[plotting_functions](https://sizespectrum.org/mizer/reference/plotting_functions.md)

Other scan functions:
[`MizerScan()`](https://sizespectrum.org/mizer/reference/MizerScan.md),
[`plotYieldVsF()`](https://sizespectrum.org/mizer/reference/plotYieldVsF.md),
[`scanEffort()`](https://sizespectrum.org/mizer/reference/scanEffort.md),
[`scanModel()`](https://sizespectrum.org/mizer/reference/scanModel.md)

## Examples

``` r
# \donttest{
scan <- scanModel(NS_params, scan_values = seq(0, 1, 0.25),
                  set_func = scanEffort(), species = c("Cod", "Herring"))
#> ℹ No `a` column so using a = 0.01 in w = a l^b, with w in g and l in cm.
#> ℹ No `b` column so using the isometric default b = 3 in w = a l^b.
plot(scan)

plot(scan, style = "envelope", mark_max = TRUE)

# }
```
