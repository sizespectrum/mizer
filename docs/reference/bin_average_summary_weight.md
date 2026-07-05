# Bin-average a summary-integral weight when second-order is enabled

Convenience wrapper around
[`bin_average_weight()`](https://sizespectrum.org/mizer/reference/bin_average_weight.md)
that is gated on the `bin_average` entry of the model's `second_order_w`
slot. When second-order bin-averaging is switched off (the default), the
weight `K` is returned unchanged so that the summary functions reproduce
the previous left-edge Riemann sums byte-for-byte. When it is switched
on, the trapezoidal bin-average of the weight is returned.

## Usage

``` r
bin_average_summary_weight(K, params)
```

## Arguments

- K:

  A numeric vector of weights indexed over the size grid, or a numeric
  matrix with the size dimension running along the columns.

- params:

  A MizerParams object whose `second_order_w` slot controls the gating.

## Value

The weight `K`, bin-averaged when
`params@second_order_w[["bin_average"]]` is `TRUE`, otherwise returned
unchanged.
