# Get default value for h

Sets `h` so that the species reaches maturity size `w_mat` at the
maturity age `age_mat` if it feeds at feeding level `f0`.

## Usage

``` r
get_h_default(params)
```

## Arguments

- params:

  A MizerParams object or a species parameter data frame

## Value

A vector with the values of h for all species

## Details

If `age_mat` is missing in the species parameter data frame, then it is
calculated from the von Bertalanffy growth curve parameters `k_vb` and
(optionally `t0`) taken from the species parameter data frame. This is
not reliable and a warning is issued.

If no growth information is given at all for a species, the default is
set to `h = 30`.

See the [Maximum Intake Rate
Coefficient](https://sizespectrum.org/mizer/articles/default_parameters.html#h-default)
section of the "Calculation of Default Parameter Values" vignette for
the mathematical derivation.

## See also

Other functions calculating defaults:
[`get_f0_default()`](https://sizespectrum.org/mizer/reference/get_f0_default.md),
[`get_gamma_default()`](https://sizespectrum.org/mizer/reference/get_gamma_default.md),
[`get_ks_default()`](https://sizespectrum.org/mizer/reference/get_ks_default.md)
