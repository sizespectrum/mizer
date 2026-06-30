# Parse the `log` argument for time-series plots

Translates the `log` argument of the time-series plotting functions into
a list of logical flags for the x and y axes, falling back to `log_x`
and `log_y` when `log` is `NULL`. For backwards compatibility a single
logical `log` value sets only `log_y`.

## Usage

``` r
parseTimePlotLog(log, log_x, log_y)
```

## Arguments

- log:

  `NULL`, a single logical value, or a character string containing only
  `"x"` and/or `"y"`.

- log_x, log_y:

  Default logical flags used when `log` is `NULL`.

## Value

A list with logical elements `log_x` and `log_y`.
