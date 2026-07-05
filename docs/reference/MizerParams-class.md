# A class to hold the parameters for a size based model.

Although it is possible to build a `MizerParams` object by hand it is
not recommended and several constructors are available. Dynamic
simulations are performed using
[`project()`](https://sizespectrum.org/mizer/reference/project.md)
function on objects of this class. As a user you should never need to
access the slots inside a `MizerParams` object directly.

## Details

The MizerParams class is fairly complex with a large number of slots,
many of which are multidimensional arrays. The dimensions of these
arrays is strictly enforced so that `MizerParams` objects are consistent
in terms of number of species and number of size classes.

The `MizerParams` class does not hold any dynamic information, e.g.
abundances or harvest effort through time. These are held in
[MizerSim](https://sizespectrum.org/mizer/reference/MizerSim-class.md)
objects.

## Slots

- `metadata`:

  A list with metadata information. See
  [`setMetadata()`](https://sizespectrum.org/mizer/reference/setMetadata.md).

- `mizer_version`:

  The package version of mizer (as returned by
  `packageVersion("mizer")`) that created or upgraded the model.

- `extensions`:

  Describes the extension chain needed to run the model. The entries are
  named by extension identifier (also the S4 marker class name) and
  ordered in S3 dispatch order, from outermost to innermost extension.
  It is either a named character vector whose values are requirement
  strings (version strings, installation specifications, or
  `NA_character_`), or a named list whose entries are length-2 character
  vectors `c(requirement = ..., version = ...)`. The `version` records
  the version of the extension package that last upgraded the object
  (`NA` if unknown) and is used by
  [`needs_upgrading()`](https://sizespectrum.org/mizer/reference/needs_upgrading.md).
  Use
  [`recordExtension()`](https://sizespectrum.org/mizer/reference/recordExtension.md)
  to write entries rather than modifying the slot directly. Extension
  subclasses are marker classes only and must not add slots.

- `time_created`:

  A POSIXct date-time object with the creation time.

- `time_modified`:

  A POSIXct date-time object with the last modified time.

- `w`:

  The size grid for the fish part of the spectrum. An increasing vector
  of weights (in grams) running from the smallest egg size to the
  largest maximum size.

- `dw`:

  The widths (in grams) of the size bins

- `w_full`:

  The size grid for the full size range including the resource spectrum.
  An increasing vector of weights (in grams) running from the smallest
  resource size to the largest maximum size of fish. The last entries of
  the vector have to be equal to the content of the w slot.

- `dw_full`:

  The width of the size bins for the full spectrum. The last entries
  have to be equal to the content of the dw slot.

- `w_min_idx`:

  A vector holding the index of the weight of the egg size of each
  species

- `maturity`:

  An array (species x size) that holds the proportion of individuals of
  each species at size that are mature. This enters in the calculation
  of the spawning stock biomass with
  [`getSSB()`](https://sizespectrum.org/mizer/reference/getSSB.md). Set
  with
  [`setReproduction()`](https://sizespectrum.org/mizer/reference/setReproduction.md).

- `psi`:

  An array (species x size) that holds the allocation to reproduction
  for each species at size, \\\psi_i(w)\\. Changed with
  [`setReproduction()`](https://sizespectrum.org/mizer/reference/setReproduction.md).

- `intake_max`:

  An array (species x size) that holds the maximum intake for each
  species at size. Changed with
  [`setMaxIntakeRate()`](https://sizespectrum.org/mizer/reference/setMaxIntakeRate.md).

- `search_vol`:

  An array (species x size) that holds the search volume for each
  species at size. Changed with
  [`setSearchVolume()`](https://sizespectrum.org/mizer/reference/setSearchVolume.md).

- `metab`:

  An array (species x size) that holds the metabolism for each species
  at size. Changed with
  [`setMetabolicRate()`](https://sizespectrum.org/mizer/reference/setMetabolicRate.md).

- `mu_b`:

  An array (species x size) that holds the external mortality rate
  \\\mu\_{ext.i}(w)\\. Changed with
  [`setExtMort()`](https://sizespectrum.org/mizer/reference/setExtMort.md).

- `ext_encounter`:

  An array (species x size) that holds the external encounter rate
  \\E\_{ext.i}(w)\\. Changed with
  [`setExtEncounter()`](https://sizespectrum.org/mizer/reference/setExtEncounter.md).

- `ext_diffusion`:

  An array (species x size) that holds the external rate at which the
  abundance density is redistributed over body size due to mixing,
  beyond the deterministic growth dynamics. Changed with
  [`ext_diffusion()`](https://sizespectrum.org/mizer/reference/setExtDiffusion.md).

- `pred_kernel`:

  An array (species x predator size x prey size) that holds the
  predation coefficient of each predator at size on each prey size. If
  this is NA then the following two slots will be used. Changed with
  [`setPredKernel()`](https://sizespectrum.org/mizer/reference/setPredKernel.md).

- `ft_pred_kernel_e`:

  An array (species x log of predator/prey size ratio) that holds the
  Fourier transform of the feeding kernel in a form appropriate for
  evaluating the encounter rate integral. If this is NA then the
  `pred_kernel` will be used to calculate the available energy integral.
  Changed with
  [`setPredKernel()`](https://sizespectrum.org/mizer/reference/setPredKernel.md).

- `ft_pred_kernel_p`:

  An array (species x log of predator/prey size ratio) that holds the
  Fourier transform of the feeding kernel in a form appropriate for
  evaluating the predation mortality integral. If this is NA then the
  `pred_kernel` will be used to calculate the integral. Changed with
  [`setPredKernel()`](https://sizespectrum.org/mizer/reference/setPredKernel.md).

- `ft_pred_kernel_d`:

  An array (species x log of predator/prey size ratio) that holds the
  Fourier transform of the feeding kernel in a form appropriate for
  evaluating the predation-diffusion integral (used when
  `use_predation_diffusion` is `TRUE`). It differs from
  `ft_pred_kernel_e` only when `second_order_w[["bin_average"]]` is
  `TRUE`, where it carries the extra power of prey size that the
  diffusion integrand needs; otherwise it equals `ft_pred_kernel_e`.
  Changed with
  [`setPredKernel()`](https://sizespectrum.org/mizer/reference/setPredKernel.md).

- `rr_pp`:

  A vector the same length as the w_full slot. The size specific growth
  rate of the resource spectrum.

- `cc_pp`:

  A vector the same length as the w_full slot. The size specific
  carrying capacity of the resource spectrum.

- `resource_dynamics`:

  Name of the function for projecting the resource abundance density by
  one timestep.

- `other_dynamics`:

  A named list of functions for projecting the values of other dynamical
  components of the ecosystem that may be modelled by a mizer extensions
  you have installed. The names of the list entries are the names of
  those components.

- `other_encounter`:

  A named list of functions for calculating the contribution to the
  encounter rate from each other dynamical component.

- `other_mort`:

  A named list of functions for calculating the contribution to the
  mortality rate from each other dynamical components.

- `other_params`:

  A list containing the parameters needed by any mizer extensions you
  may have installed to model other dynamical components of the
  ecosystem.

- `rates_funcs`:

  A named list with the names of the functions that should be used to
  calculate the rates needed by
  [`project()`](https://sizespectrum.org/mizer/reference/project.md). By
  default this will be set to the names of the built-in rate functions.

- `sc`:

  **\[experimental\]** The community abundance of the scaling community

- `species_params`:

  A data.frame to hold the species specific parameters. See
  [`species_params()`](https://sizespectrum.org/mizer/reference/species_params.md)
  for details.

- `given_species_params`:

  A data.frame to hold the species parameters that were given explicitly
  rather than obtained by default calculations.

- `gear_params`:

  Data frame with parameters for gear selectivity. See
  [`setFishing()`](https://sizespectrum.org/mizer/reference/setFishing.md)
  for details.

- `interaction`:

  The species specific interaction matrix, \\\theta\_{ij}\\. Changed
  with
  [`setInteraction()`](https://sizespectrum.org/mizer/reference/setInteraction.md).

- `selectivity`:

  An array (gear x species x w) that holds the selectivity of each gear
  for species and size, \\S\_{g,i,w}\\. Changed with
  [`setFishing()`](https://sizespectrum.org/mizer/reference/setFishing.md).

- `catchability`:

  An array (gear x species) that holds the catchability of each species
  by each gear, \\Q\_{g,i}\\. Changed with
  [`setFishing()`](https://sizespectrum.org/mizer/reference/setFishing.md).

- `initial_effort`:

  A vector containing the initial fishing effort for each gear. Changed
  with
  [`setFishing()`](https://sizespectrum.org/mizer/reference/setFishing.md).

- `initial_n`:

  An array (species x size) that holds the initial abundance of each
  species at each weight.

- `initial_n_pp`:

  A vector the same length as the w_full slot that describes the initial
  resource abundance at each weight.

- `initial_n_other`:

  A list with the initial abundances of all other ecosystem components.
  Has length zero if there are no other components.

- `resource_params`:

  List with parameters for resource.

- `A`:

  **\[deprecated\]** Formerly used to flag background species via `NA`
  values. Replaced by the `is_background` column in `species_params`.
  Will be removed in a future version.

- `linecolour`:

  A named vector of colour values, named by species. Used to give
  consistent colours in plots.

- `linetype`:

  A named vector of linetypes, named by species. Used to give consistent
  line types in plots.

- `ft_mask`:

  An array (species x w_full) with zeros for weights larger than the
  maximum weight of each species. Used to efficiently minimize
  wrap-around errors in Fourier transform calculations.

- `use_predation_diffusion`:

  A logical flag controlling whether predation-induced diffusion is
  included when calculating rates with
  [`mizerDiffusion()`](https://sizespectrum.org/mizer/reference/mizerDiffusion.md).
  Defaults to `FALSE` to preserve the behaviour of previous mizer
  versions. Set to `TRUE` to enable the diffusion term from the
  jump-growth equation.

- `second_order_w`:

  A named list with entry `flux` (the advective flux scheme: `"upwind"`,
  `"van_leer"`, or `"centred"`) and logical entry `bin_average`
  (controls whether bin-averaging is used for rates). Both default to
  the first-order setting to preserve the behaviour of previous mizer
  versions.

## See also

[`project()`](https://sizespectrum.org/mizer/reference/project.md)
[`MizerSim()`](https://sizespectrum.org/mizer/reference/MizerSim.md)
[`emptyParams()`](https://sizespectrum.org/mizer/reference/emptyParams.md)
[`newMultispeciesParams()`](https://sizespectrum.org/mizer/reference/newMultispeciesParams.md)
[`newCommunityParams()`](https://sizespectrum.org/mizer/reference/newCommunityParams.md)
[`newTraitParams()`](https://sizespectrum.org/mizer/reference/newTraitParams.md)
