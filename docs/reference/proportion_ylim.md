# Y-axis limits for a plot of a proportion

A proportion is easiest to read against the whole of the interval from 0
to 1, so that is the range a plot of one shows by default. The range is
only ever *widened* to include the data, never narrowed to the interval:
a critical feeding level or a resource level above 1 is a real feature
of the model and must stay visible.

## Usage

``` r
proportion_ylim(ylim, log_y, values)
```

## Arguments

- ylim:

  Numeric vector of length two, the limits the caller asked for.

- log_y:

  Whether the y axis is logarithmic.

- values:

  The values being plotted.

## Value

A numeric vector of length two.

## Details

Only the ends of `ylim` that the caller left as `NA` are filled in, so
an explicit limit always wins. A logarithmic axis is left alone, having
no place for the 0.
