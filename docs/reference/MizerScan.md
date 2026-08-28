# S3 class for the result of a parameter scan

**\[experimental\]**
[`scanModel()`](https://sizespectrum.org/mizer/reference/scanModel.md)
varies one aspect of a model over a range of values and measures a
quantity on the attractor the model settles on at each of them. It
returns its result as a `MizerScan` object, which is a data frame
carrying, in addition, everything that
[`plot()`](https://sizespectrum.org/mizer/reference/plot.md) needs to
draw it.

## Usage

``` r
MizerScan(
  x,
  scan_name = NULL,
  scan_units = NULL,
  value_name = NULL,
  value_units = NULL,
  type = NULL,
  params = NULL,
  reference_lines = NULL,
  settings = NULL
)

is.MizerScan(x)
```

## Arguments

- x:

  A data frame with at least three columns, laid out as described above.
  For `is.MizerScan()`, any object to test.

- scan_name:

  A string naming the quantity that was varied.

- scan_units:

  A string giving its units, for example `"1/year"`.

- value_name:

  A string naming the quantity that was measured.

- value_units:

  A string giving its units, for example `"g/year"`.

- type:

  The kind of quantity the measured values are, see
  [array_types](https://sizespectrum.org/mizer/reference/array_types.md).

- params:

  The `MizerParams` object the scan started from.

- reference_lines:

  An optional named numeric vector of x positions to mark with vertical
  lines.

- settings:

  An optional list recording the settings used.

## Value

A `MizerScan` object, which inherits from `data.frame`.

`is.MizerScan()` returns `TRUE` if `x` is a `MizerScan` object, `FALSE`
otherwise.

## Details

A `MizerScan` object behaves like an ordinary data frame with one row
for each combination of scanned value and series. Its columns are, in
order:

- 1:

  The scanned value. The column is named after the `scan_name`, so for
  example a scan over fishing effort has a column called
  `"Fishing effort"`. Use `names(scan)[[1]]` rather than hard-coding it.

- 2:

  The measured quantity, averaged over the attractor. Named after the
  `value_name`.

- `Species`:

  The series the row refers to. Named `Species` whatever the series are,
  because that is the column that mizer's colour and line-type machinery
  reads.

- `ymin`, `ymax`:

  The smallest and largest value over the sampling window. On a fixed
  point these both equal the value; on a limit cycle they give the range
  of the oscillation.

- `attractor`:

  What the state reached at this scan value is: `"fixed_point"`,
  `"limit_cycle"` or `NA` for neither. This is the column that says
  whether the value in this row can be read as an equilibrium.

- `termination`, `converged`:

  Why the run at this scan value stopped, and whether the solver met its
  own criterion. Both come from the `"convergence"` attribute that
  [`projectUntilSettled()`](https://sizespectrum.org/mizer/reference/projectUntilSettled.md)
  attaches to its result, and neither is a claim about the state — see
  `attractor` for that.

- `period`:

  The period of the limit cycle in years, or `NA`.

- `residual`:

  How far the state still is from a fixed point: the largest absolute
  relative rate of biomass change, in 1/year, taken from the
  `"convergence"` attribute that
  [`projectUntilSettled()`](https://sizespectrum.org/mizer/reference/projectUntilSettled.md)
  attaches to its result and described there. It is a biomass-weighted
  aggregate over each species' size classes, not the largest cellwise
  value of
  [`getSteadyResidual()`](https://sizespectrum.org/mizer/reference/getSteadyResidual.md).

The first three columns are the x, y and grouping variable in that
order, which is the layout
[`plotDataFrame()`](https://sizespectrum.org/mizer/reference/plotDataFrame.md)
expects.

It also carries these attributes:

- `scan_name`, `scan_units` – name and units of the quantity that was
  varied, used for the x-axis label.

- `value_name`, `value_units` – name and units of the quantity that was
  measured, used for the y-axis label.

- `type` – the kind of quantity the values are, see
  [array_types](https://sizespectrum.org/mizer/reference/array_types.md).

- `params` – the `MizerParams` object the scan started from, used for
  species colours and line types.

- `reference_lines` – an optional named numeric vector of positions on
  the x axis to mark with vertical lines, for example `c(F_MSY = 0.32)`.

- `at_max`, `max_value` – for each series, the scanned value at which
  the measured quantity is largest, and the value it takes there. See
  the section below.

- `settings` – a list recording the settings the scan was run with.

## Where the maximum is

The `at_max` attribute holds, for each series, the scanned value at
which that series' measured quantity is largest. On a
yield-versus-fishing-mortality scan that is \\F\_{MSY}\\; on a scan over
fishing effort it is the effort that maximises the quantity being
plotted. `max_value` holds the value attained there.

This is the largest value **among those that were scanned**, not the
maximum of the underlying curve. It is therefore only as good as the
grid you gave in `scan_values`, and the way to sharpen it is to scan a
finer grid near the maximum, not to interpolate a coarse one. Subsetting
a `MizerScan` with `[` recomputes both attributes from the rows that
remain, so they never go stale.

## Limitations

Because the object is a data frame subclass, its attributes survive base
R subsetting with `[` but are dropped by functions that rebuild the data
frame, including
[`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html),
[`dplyr::mutate()`](https://dplyr.tidyverse.org/reference/mutate.html),
[`subset()`](https://rdrr.io/r/base/subset.html) and
[`transform()`](https://rdrr.io/r/base/transform.html). A scan that has
lost its attributes can no longer be plotted. Subset with `[`, or
rebuild the object with `MizerScan()`.

## See also

[`scanModel()`](https://sizespectrum.org/mizer/reference/scanModel.md),
[`plot.MizerScan()`](https://sizespectrum.org/mizer/reference/plot.MizerScan.md)

Other scan functions:
[`plot.MizerScan()`](https://sizespectrum.org/mizer/reference/plot.MizerScan.md),
[`plotYieldVsF()`](https://sizespectrum.org/mizer/reference/plotYieldVsF.md),
[`scanEffort()`](https://sizespectrum.org/mizer/reference/scanEffort.md),
[`scanModel()`](https://sizespectrum.org/mizer/reference/scanModel.md)

## Examples

``` r
# \donttest{
scan <- scanModel(NS_params, scan_values = c(0, 0.5, 1),
                  set_func = scanEffort(), species = "Cod")
#> ℹ No `a` column so using a = 0.01 in w = a l^b, with w in g and l in cm.
#> ℹ No `b` column so using the isometric default b = 3 in w = a l^b.
scan
#> Biomass [g] vs Fishing effort 
#> 3 scan values x 1 series
#>  Fishing effort      Biomass Species         ymin         ymax
#>             0.0 1.218582e+12     Cod 1.218582e+12 1.218582e+12
#>             0.5 5.465367e+11     Cod 5.465367e+11 5.465367e+11
#>             1.0 3.772163e+11     Cod 3.772163e+11 3.772163e+11
#>         termination converged   attractor period     residual
#>  residual_tolerance      TRUE fixed_point     NA 0.0001509745
#>  residual_tolerance      TRUE fixed_point     NA 0.0003925230
#>  residual_tolerance      TRUE fixed_point     NA 0.0004878574
summary(scan)
#> Biomass [g] vs Fishing effort 
#> 3 scan values from 0 to 1 
#> 
#>  Species          Min          Max at_max
#>      Cod 377216256481 1.218582e+12      0
#> 
#> `at_max` is the scanned value with the largest value, over the 
#> values that were scanned. Scan a finer grid to sharpen it.
#> 
#> Attractors reached:
#> 
#> fixed_point 
#>           3 
attr(scan, "at_max")
#> Cod 
#>   0 
# }
```
