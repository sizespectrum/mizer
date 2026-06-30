# Set search volume

Set search volume

## Usage

``` r
setSearchVolume(params, search_vol = NULL, reset = FALSE, ...)

getSearchVolume(params)

search_vol(params)

search_vol(params) <- value
```

## Arguments

- params:

  MizerParams

- search_vol:

  Optional. An array (species x size) holding the search volume for each
  species at size. If not supplied, a default is set as described in the
  section "Setting search volume".

- reset:

  If set to TRUE, then the search volume will be reset to the value
  calculated from the species parameters, even if it was previously
  overwritten with a custom value. If set to FALSE (default) then a
  recalculation from the species parameters will take place only if no
  custom value has been set.

- ...:

  Unused

- value:

  search_vol

## Value

`setSearchVolume()`: A MizerParams object with updated search volume.

`getSearchVolume()` or equivalently `search_vol()`: A
`ArraySpeciesBySize` object (species x size) holding the search volume.

## Setting search volume

The search volume \\\gamma_i(w)\\ of an individual of species \\i\\ and
weight \\w\\ multiplies the predation kernel when calculating the
encounter rate in
[`getEncounter()`](https://sizespectrum.org/mizer/reference/getEncounter.md)
and the predation rate in
[`getPredRate()`](https://sizespectrum.org/mizer/reference/getPredRate.md).

The name "search volume" is a bit misleading, because \\\gamma_i(w)\\
does not have units of volume. It is simply a parameter that determines
the rate of predation. Its units depend on your choice, see section
"Units in mizer". If you have chosen to work with total abundances, then
it is a rate with units 1/year. If you have chosen to work with
abundances per m^2 then it has units of m^2/year. If you have chosen to
work with abundances per m^3 then it has units of m^3/year.

If the `search_vol` argument is not supplied, then the search volume is
set to \$\$\gamma_i(w) = \gamma_i w^q_i.\$\$ The values of \\\gamma_i\\
(the search volume at 1g) and \\q_i\\ (the allometric exponent of the
search volume) are taken from the `gamma` and `q` columns in the species
parameter dataframe. If the `gamma` column is not supplied in the
species parameter dataframe, a default is calculated by the
[`get_gamma_default()`](https://sizespectrum.org/mizer/reference/get_gamma_default.md)
function. If the `q` column is not supplied, a default of
`lambda - 2 + n` is used. Note that only for predators of size \\w = 1\\
gram is the value of the species parameter \\\gamma_i\\ the same as the
value of the search volume \\\gamma_i(w)\\.

If the `search_vol` slot has a comment and `reset = FALSE`, then a
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
[`setMaxIntakeRate()`](https://sizespectrum.org/mizer/reference/setMaxIntakeRate.md),
[`setMetabolicRate()`](https://sizespectrum.org/mizer/reference/setMetabolicRate.md),
[`setParams()`](https://sizespectrum.org/mizer/reference/setParams.md),
[`setPredKernel()`](https://sizespectrum.org/mizer/reference/setPredKernel.md),
[`setReproduction()`](https://sizespectrum.org/mizer/reference/setReproduction.md),
[`species_params()`](https://sizespectrum.org/mizer/reference/species_params.md),
[`use_predation_diffusion()`](https://sizespectrum.org/mizer/reference/use_predation_diffusion.md)

## Examples

``` r
# Inspect the current search volume
getSearchVolume(NS_params)["Cod", 1:5]
#>        0.001      0.00119      0.00142       0.0017      0.00203 
#> 6.365796e-13 7.332810e-13 8.446721e-13 9.729844e-13 1.120788e-12 

# Double the search volume for all species
sv <- getSearchVolume(NS_params) * 2
params <- setSearchVolume(NS_params, search_vol = sv)
getSearchVolume(params)["Cod", 1:5]
#>        0.001      0.00119      0.00142       0.0017      0.00203 
#> 1.273159e-12 1.466562e-12 1.689344e-12 1.945969e-12 2.241577e-12 
```
