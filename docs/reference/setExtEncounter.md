# Set external encounter rate

You will usually not need to call this function directly. Instead change
the `E_ext` and `n` species parameters with
`given_species_params(params) <-` and let mizer recalculate the external
encounter rate for you. Call `setExtEncounter()` directly only if you
want to impose a different functional form for the size dependence of
the external encounter rate. See
[`vignette("cheatsheet-changing-parameters")`](https://sizespectrum.org/mizer/articles/cheatsheet-changing-parameters.md)
for a full explanation of when to reach for which level of the model.

## Usage

``` r
setExtEncounter(params, ext_encounter = NULL, reset = FALSE, ...)

getExtEncounter(params)

ext_encounter(params)

ext_encounter(params) <- value
```

## Arguments

- params:

  MizerParams

- ext_encounter:

  Optional. An array (species x size) holding the external encounter
  rate. If not supplied, a default is calculated from the `E_ext` and
  `n` species parameters as described in the section "Setting external
  encounter rate".

- reset:

  If set to TRUE, then the external encounter rate will be reset to the
  value calculated from the species parameters, even if it was
  previously overwritten with a custom value. If set to FALSE (default)
  then a recalculation from the species parameters will take place only
  if no custom value has been set.

- ...:

  Unused

- value:

  ext_encounter

## Value

`setExtEncounter()`: A MizerParams object with updated external
encounter rate.

`getExtEncounter()` or equivalently `ext_encounter()`: A
`ArraySpeciesBySize` object (species x size) with the external encounter
rate.

## Setting external encounter rate

The external encounter rate is the rate at which a predator encounters
food that is not explicitly modelled. It is a rate with units mass/year.

The `ext_encounter` argument allows you to specify an external encounter
rate that depends on species and body size. You can see an example of
this in the Examples section of the help page for `setExtEncounter()`.

If the `ext_encounter` argument is not supplied, then the external
encounter rate is calculated as a power law: \$\$E\_{ext.i}(w) =
E\_{ext.i}\\ w^{n_i}.\$\$ The coefficient \\E\_{ext.i}\\ is taken from
the `E_ext` column of the species parameter data frame, which defaults
to 0. The exponent \\n_i\\ is taken from the `n` column of the species
parameter data frame.

If the `ext_encounter` slot has a comment and `reset = FALSE`, then a
recalculation from the species parameters is suppressed and a message is
issued if the recalculated values would differ from the stored ones.

## See also

Other functions for setting parameters:
[`gear_params()`](https://sizespectrum.org/mizer/reference/gear_params.md),
[`setExtDiffusion()`](https://sizespectrum.org/mizer/reference/setExtDiffusion.md),
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
params <- newMultispeciesParams(NS_species_params)
#> Because you have n != p, the default value for `h` is not very good.
#> Because the age at maturity is not known, I need to fall back to using
#> von Bertalanffy parameters, where available, and this is not reliable.
#> No ks column so calculating from critical feeding level.
#> Using z0 = z0pre * w_inf ^ z0exp for missing z0 values.
#> Using f0, h, lambda, kappa and the predation kernel to calculate gamma.

#### Setting allometric encounter rate #######################

# Set coefficient for each species. Here we choose 0.1 for each species
encounter_pre <- rep(0.1, nrow(species_params(params)))

# Multiply by power of size with exponent, here chosen to be 3/4
# The outer() function makes it an array species x size
allo_encounter <- outer(encounter_pre, w(params)^(3/4))

# Change the external encounter rate in the params object
ext_encounter(params) <- allo_encounter
```
