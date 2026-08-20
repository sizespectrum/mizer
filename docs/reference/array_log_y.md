# The logarithmic y axis an array's type calls for

A proportion belongs on a linear axis, so a plot of one turns `log_y`
off unless the caller asked for a particular axis. Called before
[`parsePlotLog()`](https://sizespectrum.org/mizer/reference/parsePlotLog.md),
which is why it also has to check `log`.

## Usage

``` r
array_log_y(x, log_y, log, given)
```

## Arguments

- x:

  A mizer array object.

- log_y:

  The `log_y` argument of the plot method.

- log:

  The `log` argument of the plot method.

- given:

  Whether the caller supplied `log_y` (i.e. `!missing(log_y)`).

## Value

The `log_y` to use.
