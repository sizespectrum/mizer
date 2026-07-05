# Set maximum intake rate

Set maximum intake rate

## Usage

``` r
setMaxIntakeRate(params, intake_max = NULL, reset = FALSE, ...)

getMaxIntakeRate(params)

intake_max(params)

intake_max(params) <- value
```

## Arguments

- params:

  MizerParams

- intake_max:

  Optional. An array (species x size) holding the maximum intake rate
  for each species at size. If not supplied, a default is set as
  described in the section "Setting maximum intake rate".

- reset:

  If set to TRUE, then the intake rate will be reset to the value
  calculated from the species parameters, even if it was previously
  overwritten with a custom value. If set to FALSE (default) then a
  recalculation from the species parameters will take place only if no
  custom value has been set.

- ...:

  Unused

- value:

  intake_max

## Value

`setMaxIntakeRate()`: A MizerParams object with updated maximum intake
rate.

`getMaxIntakeRate()` or equivalently `intake_max()`: A
`ArraySpeciesBySize` object (species x size) with the maximum intake
rate.

## Setting maximum intake rate

The maximum intake rate \\h_i(w)\\ of an individual of species \\i\\ and
weight \\w\\ determines the feeding level, calculated with
[`getFeedingLevel()`](https://sizespectrum.org/mizer/reference/getFeedingLevel.md).
It is measured in grams/year.

If the `intake_max` argument is not supplied, then the maximum intake
rate is set to \$\$h_i(w) = h_i w^{n_i}.\$\$ The values of \\h_i\\ (the
maximum intake rate of an individual of size 1 gram) and \\n_i\\ (the
allometric exponent for the intake rate) are taken from the `h` and `n`
columns in the species parameter dataframe. If the `h` column is not
supplied in the species parameter dataframe, it is calculated by the
[`get_h_default()`](https://sizespectrum.org/mizer/reference/get_h_default.md)
function. If the `n` column is not supplied, a default of \\n_i = 3/4\\
is used.

If \\h_i\\ is set to `Inf`, fish of species i will consume all
encountered food.

If the `intake_max` slot has a comment and `reset = FALSE`, then a
recalculation from the species parameters is suppressed and a message is
issued if the recalculated values would differ from the stored ones.

## See also

Other functions for setting parameters:
[`gear_params()`](https://sizespectrum.org/mizer/reference/gear_params.md),
[`setExtDiffusion()`](https://sizespectrum.org/mizer/reference/setExtDiffusion.md),
[`setExtEncounter()`](https://sizespectrum.org/mizer/reference/setExtEncounter.md),
[`setExtMort()`](https://sizespectrum.org/mizer/reference/setExtMort.md),
[`setFishing()`](https://sizespectrum.org/mizer/reference/setFishing.md),
[`setInteraction()`](https://sizespectrum.org/mizer/reference/setInteraction.md),
[`setMetabolicRate()`](https://sizespectrum.org/mizer/reference/setMetabolicRate.md),
[`setParams()`](https://sizespectrum.org/mizer/reference/setParams.md),
[`setPredKernel()`](https://sizespectrum.org/mizer/reference/setPredKernel.md),
[`setReproduction()`](https://sizespectrum.org/mizer/reference/setReproduction.md),
[`setSearchVolume()`](https://sizespectrum.org/mizer/reference/setSearchVolume.md),
[`species_params()`](https://sizespectrum.org/mizer/reference/species_params.md),
[`use_predation_diffusion()`](https://sizespectrum.org/mizer/reference/use_predation_diffusion.md)

## Examples

``` r
# Inspect the current maximum intake rate
getMaxIntakeRate(NS_params)["Cod", 1:5]
#>     0.001   0.00119   0.00142    0.0017   0.00203 
#> 0.6148276 0.6917271 0.7782447 0.8755836 0.9850971 

# Increase intake rate for Cod by 50%
intake_max <- getMaxIntakeRate(NS_params)
intake_max["Cod", ] <- intake_max["Cod", ] * 1.5
params <- setMaxIntakeRate(NS_params, intake_max = intake_max)
getMaxIntakeRate(params)["Cod", 1:5]
#>     0.001   0.00119   0.00142    0.0017   0.00203 
#> 0.9222414 1.0375906 1.1673671 1.3133754 1.4776457 
```
