# Calibrate the model scale to match total observed yield

**\[deprecated\]**

## Usage

``` r
calibrateYield(params, ...)
```

## Arguments

- params:

  A MizerParams object

- ...:

  Additional arguments passed to the method.

## Value

A MizerParams object. If no non-missing observed yield values are
provided, the original object is returned unchanged.

## Details

This function has been deprecated and will be removed in the future
unless you have a use case for it. If you do have a use case for it,
please let the developers know by creating an issue at
<https://github.com/sizespectrum/mizer/issues>.

Given a MizerParams object `params` for which yield observations are
available for at least some species via the `yield_observed` column in
the species_params data frame, this function returns an updated
MizerParams object which is rescaled with
[`scaleModel()`](https://sizespectrum.org/mizer/reference/scaleModel.md)
so that the total yield in the model agrees with the total observed
yield.

After using this function the total yield in the model will match the
total observed yield, summed over all species. However the yields of the
individual species will not match observations yet, with some species
having yields that are too high and others too low. So after this
function you may want to use
[`matchYields()`](https://sizespectrum.org/mizer/reference/matchYields.md).

If you have observations of species biomasses instead of yields, you can
use
[`calibrateBiomass()`](https://sizespectrum.org/mizer/reference/calibrateBiomass.md)
instead of this function.

## Examples

``` r
params <- NS_params
species_params(params)$yield_observed <-
    c(0.8, 61, 12, 35, 1.6, 20, 10, 7.6, 135, 60, 30, 78)
gear_params(params)$catchability <-
    c(1.3, 0.065, 0.31, 0.18, 0.98, 0.24, 0.37, 0.46, 0.18, 0.30, 0.27, 0.39)
params2 <- calibrateYield(params)
#> Warning: `calibrateYield()` was deprecated in mizer 2.6.0.
#> ℹ This function has not proven useful. If you do have a use case for it, please
#>   let the developers know by creating an issue at
#>   https://github.com/sizespectrum/mizer/issues
plotYieldObservedVsModel(params2)
#> The following species are not being fished in your model and will not be included in the plot: Sprat, Sandeel, N.pout.
```
