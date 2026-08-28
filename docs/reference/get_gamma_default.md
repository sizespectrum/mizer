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

The available energy is measured with the predation part of mizer's own
encounter rate,
[`mizerEncounter()`](https://sizespectrum.org/mizer/reference/mizerEncounter.md).
It is a property of the species parameters and of the resource power law
rather than of the model's dynamics, so external encounter and
contributions registered with
[`other_encounter()`](https://sizespectrum.org/mizer/reference/other_mort.md)
(including component encounter functions) are excluded. The measurement
also deliberately does not go through the extension dispatch chain, nor
through an encounter function registered with
[`setRateFunction()`](https://sizespectrum.org/mizer/reference/setRateFunction.md).

## See also

Other functions calculating defaults:
[`get_f0_default()`](https://sizespectrum.org/mizer/reference/get_f0_default.md),
[`get_h_default()`](https://sizespectrum.org/mizer/reference/get_h_default.md),
[`get_ks_default()`](https://sizespectrum.org/mizer/reference/get_ks_default.md)
