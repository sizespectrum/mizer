# Get default value for gamma

Fills in any missing values for gamma so that fish feeding on a resource
spectrum described by the power law \\\kappa w^{-\lambda}\\ achieve a
feeding level \\f_0\\. Only for internal use.

## Usage

``` r
get_gamma_default(params)
```

## Arguments

- params:

  A MizerParams object

## Value

A vector with the values of gamma for all species

## Details

See the [Search Volume
Coefficient](https://sizespectrum.org/mizer/articles/default_parameters.html#gamma-default)
section of the "Calculation of Default Parameter Values" vignette for
the mathematical derivation.

## See also

Other functions calculating defaults:
[`get_f0_default()`](https://sizespectrum.org/mizer/reference/get_f0_default.md),
[`get_h_default()`](https://sizespectrum.org/mizer/reference/get_h_default.md),
[`get_ks_default()`](https://sizespectrum.org/mizer/reference/get_ks_default.md)
