# Inspect manual rate overrides

Reports the species whose rate arrays have been set manually. The listed
parameters are retained inputs for recalculation after a reset: they no
longer determine this rate for this species, but may still contribute to
other rates or derivations. Resource overrides have `NA` in `species`.

## Usage

``` r
rate_overrides(params)
```

## Arguments

- params:

  A MizerParams object.

## Value

A data frame with one row per overridden rate and species, columns
`rate`, `species`, `parameters` (a list of parameter names), and `reset`
(a call restoring parameter-driven calculation for that species).

## See also

Other functions for setting parameters:
[`gear_params()`](https://sizespectrum.org/mizer/reference/gear_params.md),
[`setExtDiffusion()`](https://sizespectrum.org/mizer/reference/setExtDiffusion.md),
[`setExtEncounter()`](https://sizespectrum.org/mizer/reference/setExtEncounter.md),
[`setExtMort()`](https://sizespectrum.org/mizer/reference/setExtMort.md),
[`setFishing()`](https://sizespectrum.org/mizer/reference/setFishing.md),
[`setInteraction()`](https://sizespectrum.org/mizer/reference/setInteraction.md),
[`setMaxIntakeRate()`](https://sizespectrum.org/mizer/reference/setMaxIntakeRate.md),
[`setMetabolicRate()`](https://sizespectrum.org/mizer/reference/setMetabolicRate.md),
[`setParams()`](https://sizespectrum.org/mizer/reference/setParams.md),
[`setPredKernel()`](https://sizespectrum.org/mizer/reference/setPredKernel.md),
[`setReproduction()`](https://sizespectrum.org/mizer/reference/setReproduction.md),
[`setSearchVolume()`](https://sizespectrum.org/mizer/reference/setSearchVolume.md),
[`species_params()`](https://sizespectrum.org/mizer/reference/species_params.md),
[`use_predation_diffusion()`](https://sizespectrum.org/mizer/reference/use_predation_diffusion.md)

## Examples

``` r
params <- NS_params
metab(params, species = "Cod") <- 0.5 * w(params)^0.8
rate_overrides(params)
params <- setMetabolicRate(params, species = "Cod", reset = TRUE)
```
