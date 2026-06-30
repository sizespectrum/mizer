# Filter plotting data to the requested length limits

Filter plotting data to the requested length limits

## Usage

``` r
filter_plot_length_limits(plot_dat, llim)
```

## Arguments

- plot_dat:

  A data frame of plotting data, possibly with an `l` (length) column.

- llim:

  Numeric vector of length two giving the lower and upper length limits.
  Use `NA` to apply no limit at that end.

## Value

`plot_dat` filtered to the length limits, or unchanged if it has no `l`
column.
