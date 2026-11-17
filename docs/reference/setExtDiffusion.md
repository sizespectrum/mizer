# Set external diffusion rate

You will usually not need to call this function directly. Instead change
the `D_ext` and `n` species parameters with
`given_species_params(params) <-` and let mizer recalculate the external
diffusion rate for you. Call `setExtDiffusion()` directly only if you
want to impose a different functional form for the size dependence of
the external diffusion rate. See
[`vignette("guide-change-parameters")`](https://sizespectrum.org/mizer/articles/guide-change-parameters.md)
for a full explanation of when to reach for which level of the model.

## Usage

``` r
setExtDiffusion(
  params,
  ext_diffusion = NULL,
  reset = FALSE,
  ...,
  species = NULL
)

ext_diffusion(params)

ext_diffusion(params, species = NULL) <- value
```

## Arguments

- params:

  MizerParams

- ext_diffusion:

  Optional. An array (species x size) holding the external diffusion
  rate. If not supplied, a default is calculated from the `D_ext` and
  `n` species parameters as described in the section "Setting external
  diffusion rate".

- reset:

  If set to TRUE, then the external diffusion rate will be reset to the
  value calculated from the species parameters, even if it was
  previously overwritten with a custom value. If set to FALSE (default)
  then a recalculation takes place for species without a custom value.

- ...:

  Unused

- species:

  Optional character vector of species to override or reset. Omit to
  override or reset the whole array. Values may be full arrays or arrays
  containing only the selected species. For a single species its values
  without the species dimension are also accepted: a vector, or a gear x
  size matrix for the selectivity. Other species keep their existing
  override status. See
  [`rate_overrides()`](https://sizespectrum.org/mizer/reference/rate_overrides.md)
  to inspect overrides.

- value:

  ext_diffusion

## Value

`setExtDiffusion()`: A MizerParams object with updated external
diffusion rate.

`ext_diffusion()`: An `ArraySpeciesBySize` object (species x size) with
the external diffusion rate.

## Setting external diffusion rate

The external diffusion rate allows you to impose additional diffusion
beyond the predation-driven diffusion that can be internally modelled by
mizer.

The `ext_diffusion` argument allows you to specify a diffusion rate that
depends on species and body size.

If the `ext_diffusion` argument is not supplied, then the external
diffusion rate is calculated as a power law: \$\$D\_{ext.i}(w) =
D\_{ext.i}\\ w^{n_i+1}.\$\$ The coefficient \\D\_{ext.i}\\ is taken from
the `D_ext` column of the species parameter data frame, which defaults
to 0. The exponent \\n_i + 1\\ uses the `n` column of the species
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
