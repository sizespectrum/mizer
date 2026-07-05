# Get or set the second_order_w flags

**\[experimental\]** Controls whether mizer uses numerical methods that
are precise to second order in \\\Delta w\\.

## Usage

``` r
second_order_w(params)

second_order_w(params) <- value
```

## Arguments

- params:

  A MizerParams object.

- value:

  A single logical value (`TRUE` or `FALSE`) which sets both entries, a
  single flux scheme name (`"upwind"`, `"van_leer"` or `"centred"`)
  which sets only `flux`, or a named vector with entries `flux` (logical
  or scheme name) and/or `bin_average` (logical).

## Value

`second_order_w()`: A named list with entries `flux` (character) and
`bin_average` (logical).

`second_order_w<-`: A MizerParams object with the `second_order_w` flags
updated and, when `bin_average` is changed, all model parameters
recalculated via
[`setParams()`](https://sizespectrum.org/mizer/reference/setParams.md).

## Details

The slot is a named list with entries:

- `flux`:

  The advective-flux reconstruction scheme used in the numerical solver.
  `"upwind"` is the first-order upwind scheme. `"van_leer"` is the
  second-order scheme with the total-variation- diminishing van Leer
  limiter, which keeps abundances non-negative. `"centred"` is the
  second-order scheme with the unlimited centred flux, which is
  genuinely second order even at extrema but is not
  monotonicity-preserving (it can produce small over/undershoots and is
  best used with some physical diffusion).

- `bin_average`:

  Logical. Controls whether bin-averaging is used for quantities that
  need it in order to be second-order precise in bin size. When `FALSE`,
  point-sampling at the left bin edge is used.

When `flux` is `"upwind"` and `bin_average` is `FALSE` (the defaults),
mizer preserves the behaviour of previous mizer versions. Setting both
to their second-order values gives a consistently second-order model.

The setter accepts a single logical value (which sets both entries), a
single scheme name (which sets only `flux`), or a named vector to set
individual entries. The setter re-runs
[`setParams()`](https://sizespectrum.org/mizer/reference/setParams.md)
to rebuild precomputed arrays when `bin_average` is changed.
