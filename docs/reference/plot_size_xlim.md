# Choose the x-axis limits for a given size axis

Choose the x-axis limits for a given size axis

## Usage

``` r
plot_size_xlim(wlim, size_axis, llim = c(NA, NA))
```

## Arguments

- wlim:

  Numeric vector of length two giving the weight limits.

- size_axis:

  Either `"w"` (weight) or `"l"` (length).

- llim:

  Numeric vector of length two giving the length limits.

## Value

`llim` for a length axis, otherwise `wlim`.
