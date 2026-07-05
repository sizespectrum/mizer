# Helper function to produce nice breaks on logarithmic axes

This is needed when the logarithmic y-axis spans less than one order of
magnitude, in which case the ggplot2 default produces no ticks.

## Usage

``` r
log_breaks(n = 6)
```

## Arguments

- n:

  Approximate number of ticks

## Value

A function that can be used as the break argument in calls to
scale_y_continuous() or scale_x_continuous()

## Details

Thanks to Heather Turner at
<https://stackoverflow.com/questions/14255533/pretty-ticks-for-log-normal-scale-using-ggplot2-dynamic-not-manual>
