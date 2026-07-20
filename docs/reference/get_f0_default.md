# Get default value for f0

Fills in any missing values for f0 so that if the prey abundance was
described by the power law \\\kappa w^{-\lambda}\\ then the encounter
rate coming from the given `gamma` parameter would lead to the feeding
level \\f_0\\. This is thus doing the inverse of
[`get_gamma_default()`](https://sizespectrum.org/mizer/reference/get_gamma_default.md).
Only for internal use.

## Usage

``` r
get_f0_default(params)
```

## Arguments

- params:

  A MizerParams object

## Value

A vector with the values of f0 for all species

## Details

For species for which no value for `gamma` is specified in the species
parameter data frame, the `f0` values is kept as provided in the species
parameter data frame or it is set to 0.6 if it is not provided.

See the [Target Feeding
Level](https://sizespectrum.org/mizer/articles/default_parameters.html#f0-default)
section of the "Calculation of Default Parameter Values" vignette for
the mathematical derivation.

## See also

Other functions calculating defaults:
[`get_gamma_default()`](https://sizespectrum.org/mizer/reference/get_gamma_default.md),
[`get_h_default()`](https://sizespectrum.org/mizer/reference/get_h_default.md),
[`get_ks_default()`](https://sizespectrum.org/mizer/reference/get_ks_default.md)
