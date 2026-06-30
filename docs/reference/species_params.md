# Species parameters

These functions allow you to get or set the species-specific parameters
stored in a MizerParams object.

## Usage

``` r
species_params(params)

species_params(params) <- value

given_species_params(params)

given_species_params(params) <- value

calculated_species_params(params)
```

## Arguments

- params:

  A MizerParams object

- value:

  A data frame with the species parameters

## Value

`species_params()`: Data frame containing all species parameters
currently stored in the model.

`species_params<-()`: Updates the full species parameter table after
validating it with
[`validSpeciesParams()`](https://sizespectrum.org/mizer/reference/validSpeciesParams.md)
and then recalculating the model parameters with
[`setParams()`](https://sizespectrum.org/mizer/reference/setParams.md).

`given_species_params()`: Data frame containing the species parameter
values that were supplied explicitly by the user.

`given_species_params<-()`: Updates the explicitly supplied species
parameters after validating them with
[`validGivenSpeciesParams()`](https://sizespectrum.org/mizer/reference/validSpeciesParams.md)
and then recalculating the full species parameter table and dependent
model quantities.

`calculated_species_params()`: Data frame containing only those species
parameter entries that are not explicit user input. Columns that would
consist entirely of `NA` values are dropped.

## Details

There are a lot of species parameters and we will list them all below,
but most of them have sensible default values. The only required columns
are `species` for the species name and `w_inf` for its von Bertalanffy
asymptotic size. However if you have information about the values of
other parameters then you should provide them.

Three species parameters describe maximum sizes and play distinct roles:

- `w_inf` is the von Bertalanffy asymptotic size of an average
  individual. It is the required maximum-size parameter and is used to
  set default values for `w_max`, `w_repro_max` and `w_mat`.

- `w_repro_max` is the size at which a typical mature individual invests
  all of its available energy into reproduction, see
  [`setReproduction()`](https://sizespectrum.org/mizer/reference/setReproduction.md).
  It is not a hard ceiling on size and defaults to `w_inf`.

- `w_max` is purely a computational boundary: it sets the upper end of
  the size grid and the range of plots. It defaults to `1.5 * w_inf`.
  For backwards compatibility, if `w_inf` is not supplied it is taken
  from `w_repro_max` or `w_max` instead.

Mizer distinguishes between the species parameters that you have given
explicitly and the species parameters that have been calculated by mizer
or set to default values. You can retrieve the given species parameters
with `given_species_params()` and the calculated ones with
`calculated_species_params()`. You get all species_params with
`species_params()`.

If you change given species parameters with `given_species_params<-()`
this will trigger a re-calculation of the calculated species parameters,
where necessary. However if you change species parameters with
`species_params<-()` no recalculation will take place and furthermore
your values could be overwritten by a future recalculation triggered by
a call to `given_species_params<-()` . So in most use cases you will
only want to use `given_species_params<-()`.

There are some species parameters that are used to set up the
size-dependent parameters that are used in the mizer model:

- `gamma` and `q` are used to set the search volume, see
  [`setSearchVolume()`](https://sizespectrum.org/mizer/reference/setSearchVolume.md).

- `h` and `n` are used to set the maximum intake rate, see
  [`setMaxIntakeRate()`](https://sizespectrum.org/mizer/reference/setMaxIntakeRate.md).

- `k`, `ks` and `p` are used to set activity and basic metabolic rate,
  see
  [`setMetabolicRate()`](https://sizespectrum.org/mizer/reference/setMetabolicRate.md).

- `z0`, `z_ext` and `d` are used to set the external mortality rate, see
  [`setExtMort()`](https://sizespectrum.org/mizer/reference/setExtMort.md).

- `E_ext` and `n` are used to set the external encounter rate, see
  [`setExtEncounter()`](https://sizespectrum.org/mizer/reference/setExtEncounter.md).

- `D_ext` and `n` are used to set the external diffusion rate, see
  [`setExtDiffusion()`](https://sizespectrum.org/mizer/reference/setExtDiffusion.md).

- `w_mat`, `w_mat25`, `w_repro_max` and `m` are used to set the
  allocation to reproduction, see
  [`setReproduction()`](https://sizespectrum.org/mizer/reference/setReproduction.md).

- `pred_kernel_type` specifies the shape of the predation kernel. The
  default is a "lognormal", for other options see the "Setting predation
  kernel" section in the help for
  [`setPredKernel()`](https://sizespectrum.org/mizer/reference/setPredKernel.md).

- `beta` and `sigma` are parameters of the lognormal predation kernel,
  see
  [`lognormal_pred_kernel()`](https://sizespectrum.org/mizer/reference/lognormal_pred_kernel.md).
  There will be other parameters if you are using other predation kernel
  functions.

When you change one of the above species parameters using
`given_species_params<-()` or `species_params<-()`, the new value will
be used to update the corresponding size-dependent rates automatically,
unless you have set those size-dependent rates manually, in which case
the corresponding species parameters will be ignored.

There are some species parameters that are used directly in the model
rather than being used for setting up size-dependent parameters:

- `alpha` is the assimilation efficiency, the proportion of the consumed
  biomass that can be used for growth, metabolism and reproduction, see
  the help for
  [`getEReproAndGrowth()`](https://sizespectrum.org/mizer/reference/getEReproAndGrowth.md).

- `w_min` is the egg size.

- `interaction_resource` sets the interaction strength with the
  resource, see "Predation encounter" section in the help for
  [`getEncounter()`](https://sizespectrum.org/mizer/reference/getEncounter.md).

- `erepro` is the reproductive efficiency, the proportion of the energy
  invested into reproduction that is converted to egg biomass, see
  [`getRDI()`](https://sizespectrum.org/mizer/reference/getRDI.md).

- `R_max` is the parameter in the Beverton-Holt density dependence added
  to the reproduction, see
  [`setBevertonHolt()`](https://sizespectrum.org/mizer/reference/setBevertonHolt.md).
  There will be other such parameters if you use other density
  dependence functions, see the "Density dependence" section in the help
  for
  [`setReproduction()`](https://sizespectrum.org/mizer/reference/setReproduction.md).

Two parameters are used only by functions that need to convert between
weight and length:

- `a` and `b` are the parameters in the allometric weight-length
  relationship \\w = a l ^ b\\.

If you have supplied the `a` and `b` parameters, then you can replace
weight parameters like `w_inf`, `w_max`, `w_mat`, `w_mat25`,
`w_repro_max` and `w_min` by their corresponding length parameters
`l_inf`, `l_max`, `l_mat`, `l_mat25`, `l_repro_max` and `l_min`.

The parameters that are only used to calculate default values for other
parameters are:

- `f0` is the feeding level and is used to get a default value for the
  coefficient of the search volume `gamma`, see
  [`get_gamma_default()`](https://sizespectrum.org/mizer/reference/get_gamma_default.md).

- `fc` is the critical feeding level below which the species can not
  maintain itself. This is used to get a default value for the
  coefficient `ks` of the metabolic rate, see
  [`get_ks_default()`](https://sizespectrum.org/mizer/reference/get_ks_default.md).

- `age_mat` is the age at maturity and is used to get a default value
  for the coefficient `h` of the maximum intake rate, see
  [`get_h_default()`](https://sizespectrum.org/mizer/reference/get_h_default.md).

- If `age_mat` is not supplied, mizer used the von Bertalanffy
  parameters `k_vb`, `w_inf` and `t0` as well as the weight-length
  exponent `b` to determine it. This is unreliable and is therefore not
  recommended.

Changing these parameters with `species_params<-()` updates the stored
species parameter table and triggers a recalculation via
[`setParams()`](https://sizespectrum.org/mizer/reference/setParams.md).
However they only affect model behaviour if the corresponding downstream
parameters are recalculated rather than kept at explicitly supplied
values. In typical workflows these quantities should therefore be
changed via `given_species_params<-()`.

There are other species parameters that are used in tuning the model to
observations:

- `biomass_observed` and `biomass_cutoff` allow you to specify for each
  species the total observed biomass above some cutoff size. This is
  used by
  [`calibrateBiomass()`](https://sizespectrum.org/mizer/reference/calibrateBiomass.md)
  and
  [`matchBiomasses()`](https://sizespectrum.org/mizer/reference/matchBiomasses.md).

- `yield_observed` allows you to specify for each species the total
  annual fisheries yield. This is used by
  [`calibrateYield()`](https://sizespectrum.org/mizer/reference/calibrateYield.md)
  and
  [`matchYields()`](https://sizespectrum.org/mizer/reference/matchYields.md).

Finally there are two species parameters that control the way the
species are represented in plots:

- `linecolour` specifies the colour and can be any valid R colour value.

- `linetype` specifies the line type ("solid", "dashed", "dotted",
  "dotdash", "longdash", "twodash" or "blank")

Other species-specific information that is related to how the species is
fished is specified in a gear parameter data frame, see
[`gear_params()`](https://sizespectrum.org/mizer/reference/gear_params.md).
However in the case where each species is caught by only a single gear,
this information can also optionally be provided as species parameters
and
[`newMultispeciesParams()`](https://sizespectrum.org/mizer/reference/newMultispeciesParams.md)
will transfer them to the `gear_params` data frame. However changing
these parameters later in the species parameter data frames will have no
effect.

You are allowed to include additional columns in the species parameter
data frames. They will simply be ignored by mizer but will be stored in
the MizerParams object, in case your own code makes use of them.

## See also

[`validSpeciesParams()`](https://sizespectrum.org/mizer/reference/validSpeciesParams.md),
[`setParams()`](https://sizespectrum.org/mizer/reference/setParams.md)

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
[`use_predation_diffusion()`](https://sizespectrum.org/mizer/reference/use_predation_diffusion.md)
