# Set external encounter rate

You will usually not need to call this function directly. Instead change
the `E_ext` and `n` species parameters with
`given_species_params(params) <-` and let mizer recalculate the external
encounter rate for you. Call `setExtEncounter()` directly only if you
want to impose a different functional form for the size dependence of
the external encounter rate. See
[`vignette("guide-change-parameters")`](https://sizespectrum.org/mizer/articles/guide-change-parameters.md)
for a full explanation of when to reach for which level of the model.

## Usage

``` r
setExtEncounter(
  params,
  ext_encounter = NULL,
  reset = FALSE,
  ...,
  species = NULL
)

ext_encounter(params)

ext_encounter(params, species = NULL) <- value
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
  then a recalculation takes place for species without a custom value.

- ...:

  Unused

- species:

  Optional character vector of species to override or reset. Omit to
  override or reset the whole array. Values may be full arrays or arrays
  containing only the selected species; for one species a vector is also
  accepted. Other species keep their existing override status. See
  [`rate_overrides()`](https://sizespectrum.org/mizer/reference/rate_overrides.md)
  to inspect overrides.

- value:

  ext_encounter

## Value

`setExtEncounter()`: A MizerParams object with updated external
encounter rate.

`ext_encounter()`: An `ArraySpeciesBySize` object (species x size) with
the external encounter rate.

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

With `reset = FALSE`, manually overridden species are preserved while
the other species are recalculated. A message is issued if recalculation
would change protected values. See
[`rate_overrides()`](https://sizespectrum.org/mizer/reference/rate_overrides.md)
to inspect overrides.

## See also

Other functions for setting parameters:
[`gear_params()`](https://sizespectrum.org/mizer/reference/gear_params.md),
[`rate_overrides()`](https://sizespectrum.org/mizer/reference/rate_overrides.md),
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

#### Setting allometric encounter rate #######################

# Set coefficient for each species. Here we choose 0.1 for each species
encounter_pre <- rep(0.1, nrow(species_params(params)))

# Multiply by power of size with exponent, here chosen to be 3/4
# The outer() function makes it an array species x size
allo_encounter <- outer(encounter_pre, w(params)^(3/4))

# Change the external encounter rate in the params object
ext_encounter(params) <- allo_encounter
```
