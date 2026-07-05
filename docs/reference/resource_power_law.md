# Construct the background resource power-law spectrum

Internal helper returning the auto-calculated resource power law
\\\kappa\\ w^{-\lambda}\\ on the full size grid, optionally truncated at
an upper cutoff `w_max`. When the `bin_average` entry of the model's
`second_order_w` slot is set, the exact bin average of the power law
over each bin is returned instead (with the bin straddling `w_max`
getting the partial average), so that the initial resource is the
finite-volume cell average of the background spectrum, consistent with
the bin-integrated encounter convolution that consumes it as a cell
average. Otherwise the left-edge point values are returned,
byte-identical to previous mizer.

## Usage

``` r
resource_power_law(params, kappa, lambda, w_max = Inf)
```

## Arguments

- params:

  A MizerParams object whose `second_order_w` slot controls the gating
  and whose `w_full`/`dw_full` give the grid.

- kappa:

  The coefficient \\\kappa\\ of the power law.

- lambda:

  The exponent so the power law is \\w^{-\lambda}\\.

- w_max:

  Optional upper cutoff. The power law is taken to be zero at and above
  `w_max`. Defaults to `Inf` (no cutoff).

## Value

A numeric vector (same length as `w_full`) of the resource spectrum.

## Details

This is used wherever the background resource spectrum \\\kappa
w^{-\lambda}\\ is constructed from scratch: the initial resource
abundance in
[`newMultispeciesParams()`](https://sizespectrum.org/mizer/reference/newMultispeciesParams.md)
and the temporary prey spectra used to compute the default `gamma`/`f0`
([`get_gamma_default()`](https://sizespectrum.org/mizer/reference/get_gamma_default.md),
[`get_f0_default()`](https://sizespectrum.org/mizer/reference/get_f0_default.md))
and the consumer initial abundances
([`get_initial_n()`](https://sizespectrum.org/mizer/reference/get_initial_n.md)).
The bin-averaged resource capacity and rate are produced directly by
[`setResource()`](https://sizespectrum.org/mizer/reference/setResource.md).
