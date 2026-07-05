# Get default value for `ks`

Fills in any missing values for `ks` so that the critical feeding level
needed to sustain the species is as specified in the `fc` column in the
species parameter data frame. If that column is not provided the default
critical feeding level \\f_c = 0.2\\ is used.

## Usage

``` r
get_ks_default(params)
```

## Arguments

- params:

  A MizerParams object

## Value

A vector with the values of ks for all species

## See also

Other functions calculating defaults:
[`get_f0_default()`](https://sizespectrum.org/mizer/reference/get_f0_default.md),
[`get_gamma_default()`](https://sizespectrum.org/mizer/reference/get_gamma_default.md),
[`get_h_default()`](https://sizespectrum.org/mizer/reference/get_h_default.md)
