# Species parameters

These functions allow you to get or set the species-specific parameters
stored in a MizerParams object.

## Usage

``` r
species_params(object, ...)

species_params(object, recalculate = TRUE) <- value

is.species_params(x)

given_species_params(object, ...)

is.given_species_params(x)

given_species_params(object) <- value

calculated_species_params(params)
```

## Arguments

- object:

  A MizerParams object, a MizerSim object or a data frame

- ...:

  Other arguments passed to methods.

- strict:

  Whether to raise an error, rather than correct silently, for the
  inconsistencies that can be corrected. Used internally.

- check_misspellings:

  Whether to report column names that look like misspellings of standard
  species parameter names. `TRUE` by default. Mizer passes `FALSE` when
  it re-validates a table whose columns it has already checked, so that
  the report is made once, when the column is introduced, rather than
  again on every rebuild.

- recalculate:

  Whether `species_params<-()` should be allowed to re-derive calculated
  species parameters and rates that depend on a changed parameter.
  Defaults to `TRUE`; mizer still skips the rebuild when all changes are
  to columns with no cached dependants. See the section "Setting species
  parameters without recalculation" below before setting it to `FALSE`.

- value:

  A data frame with the new species parameters.

- x:

  An object to test with `is.species_params()` or
  `is.given_species_params()`.

- params:

  A MizerParams object.

## Value

`species_params()`: Data frame containing all species parameters
currently stored in the model.

`species_params<-()`: Updates the `given_species_params` with any
parameters you have changed, and recalculates the full species parameter
table and model parameters when a changed column has cached dependants.
A column that `value` does not have is one you no longer supply and is
removed from the given species parameters, see the section "Removing a
species parameter column" below. With `recalculate = FALSE` it only does
the recording and stores the parameters you supplied, see the section
"Setting species parameters without recalculation" below.

`given_species_params()`: Data frame containing the species parameter
values that were supplied explicitly by the user.

`given_species_params<-()`: Replaces the authoritative table of
parameters that are to count as explicit user input. Every non-`NA`
entry in `value` is recorded as given, even when it is numerically equal
to the value currently in `species_params()`. This lets you protect a
calculated value against future recalculation. An `NA` entry, or removal
of a column, hands a previously given parameter back to mizer's
calculation – and where there is no mizer calculation to hand it back
to, removing the column removes the parameter from the model, see the
section "Removing a species parameter column" below. Dependent
quantities are recalculated only when the replacement can change them;
merely marking the current value as given does not rebuild the model.

This setter also warns when a change you asked for cannot take effect,
namely when the parameter is overridden by another one you have already
given (`f0` by `gamma`, `fc` by `ks`, `age_mat` by `h`), when the rate
array it feeds has been set by hand and so is no longer calculated, or
when it is a gear parameter that mizer reads from
[`gear_params()`](https://sizespectrum.org/mizer/reference/gear_params.md)
instead. `species_params<-()` stays quiet about all three.
`given_species_params<-()` has no `recalculate` argument; where you need
to record values without recalculating, use `species_params<-()` or
[`record_given_species_params()`](https://sizespectrum.org/mizer/reference/record_given_species_params.md).

`calculated_species_params()`: Data frame containing only those species
parameter entries that are not explicit user input. Columns that would
consist entirely of `NA` values are dropped.

`is.species_params()` returns `TRUE` if `x` is a `species_params`
object, `FALSE` otherwise.

`is.given_species_params()` returns `TRUE` if `x` is a
`given_species_params` object, `FALSE` otherwise.

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

When you change species parameters with `species_params<-()`, mizer
automatically detects which parameters you have changed. It records
these changed parameters in `given_species_params` so that they are
protected against being overwritten by future recalculations. It then
re-calculates the quantities that depend on the changed parameters.
Changes to observation, direct-runtime or other custom columns that base
mizer does not use to build a cached quantity do not trigger that
recalculation. Unknown columns on an extension object retain the
conservative recalculation path because an extension setter may use
them.

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
  The Gaussian-mixture kernel instead uses the list-columns `kernel_p`,
  `kernel_mean`, and `kernel_sd`, see
  [`gaussian_mixture_pred_kernel()`](https://sizespectrum.org/mizer/reference/gaussian_mixture_pred_kernel.md).
  There will be other parameters if you are using other predation kernel
  functions.

When you change one of the above species parameters using
`species_params<-()` or `given_species_params<-()`, the new value will
be used to update the corresponding size-dependent rates automatically,
unless you have set those size-dependent rates manually, in which case
the corresponding species parameters will be ignored. Mizer warns you
when that happens, because the value is then in the species parameter
table without having any effect on the model. The warning names the call
that puts the rate back under the control of the species parameters.

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

You can also keep both, and change either of them later. Mizer keeps the
two consistent by the rule that the one you gave last wins, and if you
gave both at the same time the weight wins. So on a model set up with
lengths you can still set `w_mat` with `species_params<-()` and mizer
will update `l_mat` to match, and if you set `l_mat` it will update
`w_mat` as always. When you supply a length and a weight together that
do not agree, mizer uses the weight and warns you that it has changed
the length to match.

The rule is applied when a species parameter data frame is put into a
model. A data frame that you have taken out of a model and are editing
on its own is left exactly as you write it: no conversions, no checks
and no warnings until you assign it back with `species_params<-()` or
`given_species_params<-()`, which is when mizer can tell which values
you changed. A data frame that was never in a model, for example one you
pass to
[`validSpeciesParams()`](https://sizespectrum.org/mizer/reference/validSpeciesParams.md),
carries no such history, so a length and a weight that disagree there
count as given at the same time and the weight wins.

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

Changing these parameters with `species_params<-()` will trigger a
recalculation of the downstream parameters, provided they are not
protected by being explicitly given.

There are other species parameters that are used in tuning the model to
observations:

- `biomass_observed` and `biomass_cutoff` allow you to specify for each
  species the total observed biomass above some cutoff size. This is
  used by
  [`calibrateBiomass()`](https://sizespectrum.org/mizer/reference/calibrateBiomass.md)
  and
  [`matchBiomasses()`](https://sizespectrum.org/mizer/reference/matchBiomasses.md).

The total annual fisheries yield is not a species parameter but a gear
parameter, because it is observed for each gear separately, see
[`gear_params()`](https://sizespectrum.org/mizer/reference/gear_params.md).
For backwards compatibility mizer still accepts a `yield_observed`
column in the species parameter data frame, see
[`get_yield_observed()`](https://sizespectrum.org/mizer/reference/get_yield_observed.md).

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

## Extracting a column with `$`

`species_params(params)$w_mat` returns the column as a vector named
after the species. Unlike `$` on an ordinary data frame, it does **not**
partially match the column name. Partial matching is dangerous here
because so many species parameter names are prefixes of others: in a
model without length-weight parameters `species_params(params)$a` used
to return the `alpha` column and `$b` the `beta` column, complete with
species names, so code converting weights to lengths silently got the
assimilation efficiency and the preferred predator/prey mass ratio
instead. Writing was never partially matched, so reads and writes
disagreed about which column `$b` meant.

A name that is not a column now gives `NULL`, so
`is.null(species_params(params)$foo)` is a reliable way of testing
whether a parameter is present. If the name would have partially matched
a column under the old behaviour you also get a warning naming that
column, because that is exactly the case where existing code changes its
meaning. The same holds for
[`gear_params()`](https://sizespectrum.org/mizer/reference/gear_params.md).

## Removing a species parameter column

A column that is missing from the table you assign is one you no longer
supply, and mizer takes it out of the given species parameters. What
happens next depends on whether mizer knows how to produce the parameter
itself: one that it calculates comes straight back as a calculated
value, while one that it knows nothing about leaves the model
altogether. Both setters work this way, so either of these removes a
custom column:

    species_params(params)$my_col <- NULL
    given_species_params(params)$my_col <- NULL

and either of these hands a parameter you had given back to mizer,
exactly as setting it to `NA` in the given species parameters does:

    species_params(params)$gamma <- NULL
    given_species_params(params)$gamma <- NULL

This is how an extension package withdraws a species parameter it added
when the user switches the extension off, without reaching into the
`species_params` slot. Mizer reports the removal at `info_level` 3, see
[`setParams()`](https://sizespectrum.org/mizer/reference/setParams.md).

Because the whole table is compared, assigning a table with only some of
the model's columns withdraws all the others. There is not much room for
surprise: `species_params<-()` validates what you give it, so a table
without `species` and one of `w_inf`, `w_max` or `w_repro_max` is an
error rather than a partial update. Still, edit the table you got from
the accessor rather than building a new one from a handful of columns.

## Setting species parameters without recalculation

`species_params(params, recalculate = FALSE) <- value` records the
values you changed among the given species parameters, so that they are
not calculated away later, and stores `value` as the species parameters.
It then stops there: the calculated species parameters are not
re-derived from the given ones, no missing parameters are filled in with
their default values, and none of the size-dependent rates are
recalculated. Your species parameters are stored as you supplied them,
after the same checks and length-to-weight conversions that writing into
the `species_params` slot would trigger.

This is for code that has worked out a species parameter *together with*
the rate array that the parameter determines, for example an optimiser
that fits `ks` and the matching `metab`, or `z_ext` and the matching
`mu_b`. There the recalculation is not just wasted work but would
overwrite the rates the caller has just set.

The object you get back is only as consistent as you make it. Mizer will
not check that the species parameters you supplied agree with the rate
arrays in the object, nor that they agree with the other species
parameters that are normally derived from them. Unless you are setting
the affected rates yourself, use the default `recalculate = TRUE`.

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
