# Set metabolic rate

Sets the rate at which energy is used for metabolism and activity

## Usage

``` r
setMetabolicRate(object, metab = NULL, p = NULL, reset = FALSE, ...)

getMetabolicRate(params)

metab(params)

metab(params) <- value
```

## Arguments

- object:

  A MizerParams object

- metab:

  Optional. An array (species x size) holding the metabolic rate for
  each species at size. If not supplied, a default is set as described
  in the section "Setting metabolic rate".

- p:

  The allometric metabolic exponent. This is only used if `metab` is not
  given explicitly and if the exponent is not specified in a `p` column
  in the `species_params`.

- reset:

  If set to TRUE, then the metabolic rate will be reset to the value
  calculated from the species parameters, even if it was previously
  overwritten with a custom value. If set to FALSE (default) then a
  recalculation from the species parameters will take place only if no
  custom value has been set.

- ...:

  Unused

- params:

  A MizerParams object

- value:

  metab

## Value

`setMetabolicRate()`: A MizerParams object with updated metabolic rate.

`getMetabolicRate()` or equivalently `metab()`: A `ArraySpeciesBySize`
object (species x size) with the metabolic rate.

## Setting metabolic rate

The metabolic rate is subtracted from the energy income rate to
calculate the rate at which energy is available for growth and
reproduction, see
[`getEReproAndGrowth()`](https://sizespectrum.org/mizer/reference/getEReproAndGrowth.md).
It is measured in grams/year.

If the `metab` argument is not supplied, then for each species the
metabolic rate \\k(w)\\ for an individual of size \\w\\ is set to
\$\$k(w) = k_s w^p + k w,\$\$ where \\k_s w^p\\ represents the rate of
standard metabolism and \\k w\\ is the rate at which energy is expended
on activity and movement. The values of \\k_s\\, \\p\\ and \\k\\ are
taken from the `ks`, `p` and `k` columns in the species parameter
dataframe. If any of these parameters are not supplied, the defaults are
\\k = 0\\, \\p = 3/4\\ and \$\$k_s = f_c h \alpha w\_{mat}^{n-p},\$\$
where \\f_c\\ is the critical feeding level taken from the `fc` column
in the species parameter data frame. If the critical feeding level is
not specified, a default of \\f_c = 0.2\\ is used.

If the `metab` slot has a comment and `reset = FALSE`, then a
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
[`setParams()`](https://sizespectrum.org/mizer/reference/setParams.md),
[`setPredKernel()`](https://sizespectrum.org/mizer/reference/setPredKernel.md),
[`setReproduction()`](https://sizespectrum.org/mizer/reference/setReproduction.md),
[`setSearchVolume()`](https://sizespectrum.org/mizer/reference/setSearchVolume.md),
[`species_params()`](https://sizespectrum.org/mizer/reference/species_params.md),
[`use_predation_diffusion()`](https://sizespectrum.org/mizer/reference/use_predation_diffusion.md)

## Examples

``` r
# Inspect the current metabolic rate
getMetabolicRate(NS_params)["Cod", 1:5]
#>      0.001    0.00119    0.00142     0.0017    0.00203 
#> 0.09767498 0.11054111 0.12510202 0.14158096 0.16023056 

# Reset metabolic rate from species parameters
params <- setMetabolicRate(NS_params, reset = TRUE)
getMetabolicRate(params)["Cod", 1:5]
#>      0.001    0.00119    0.00142     0.0017    0.00203 
#> 0.09767498 0.11054111 0.12510202 0.14158096 0.16023056 
```
