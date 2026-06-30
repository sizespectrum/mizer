# Parse the log-axis arguments of a mizer plot function

Internal helper that resolves the various ways of specifying which axes
should use a logarithmic scale into a consistent pair of logical flags.
It is exported so that extension packages (such as mizerMR) can reuse it
in their own array
[`plot()`](https://sizespectrum.org/mizer/reference/plot.md) methods.

## Usage

``` r
parsePlotLog(log, log_x = FALSE, log_y = FALSE)
```

## Arguments

- log:

  Either `NULL`, a single logical (legacy form, toggling only the
  y-axis), or a character string containing only the letters `"x"`
  and/or `"y"` to indicate which axes should be logarithmic.

- log_x, log_y:

  Default logical flags used when `log` is `NULL`.

## Value

A list with logical components `log_x` and `log_y`.
