# Restrict plot data to a range of weights

Internal helper that filters a plot data frame to the weight range given
by `wlim`. It is exported so that extension packages (such as mizerMR)
can reuse it in their own array
[`plot()`](https://sizespectrum.org/mizer/reference/plot.md) methods.

## Usage

``` r
apply_wlim(data, wlim)
```

## Arguments

- data:

  A data frame with a numeric `w` column.

- wlim:

  A length-2 numeric vector giving the lower and upper weight limits.
  Either entry may be `NA` to leave that side unrestricted.

## Value

The subset of `data` with `w` inside `wlim`.
