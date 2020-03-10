
# mizer 1.0.1.9000 

## Backwards compatibility

This version of mizer introduces a large number of new features but is almost
fully backwards compatible with version 1.0 with the exception of bug fixes and
the following breaking changes:

* Removed the `print_it` argument from plot functions.
* plotFeedingLevel() now only plots the values within the size range of each
  species. If for some reason you want the old plots that show a feeding level
  also for sizes that the fish can never have, you need to supply an argument
  `all.sizes = TRUE`.
* The way the density-dependence in the reproduction rate is set has changed,
  see `RDD` argument in `setReproduction()`.
* The `sex_ratio` argument has been removed from `getRDI()` and `getRDD()`.
* The `set_scaling_model()` function has been removed because such models can
  now be set up with `newTraitParams()` with the options `perfect_scaling = TRUE`
  and `egg_size_scaling = TRUE`.
* The functions `display_frames()`, `addSpecies()`, `setBackground()` and 
  `retuneAbundance()` have been removed to the "mizerExperimental" package
  (https://sizespectrum.org/mizerExperimental)
* During runs of `project()` a progress bar is displayed by default. You can 
  turn this off with the option `progress_bar = FALSE.

## Setting up models

The new functions

* `newCommunityParams()`
* `newTraitParams()`
* `newMultispeciesParams()`

replace the old functions `set_community_model()`, `set_trait_model()` and
`MizerParams()`, which are now deprecated. The new functions choose better
default values, in particular for metabolic rate and maximum intake rate.

## Setting model parameters
After setting up a mizer model, it is possible to change specific model
parameters with the new functions

* `setPredationKernel()`
* `setSearchVolume()`
* `setInteraction()`
* `setMaxIntakeRate()`
* `setMetabolicRate()`
* `setExtMortality()`
* `setReproduction()`
* `setFishing()`
* `setPlankton()`

The new function `setParams()` is a wrapper for all of the above functions
and is also used when setting up a new model with `newMultispeciesParams()`.
(#51)

The documentation for these functions serves to explain the details of the
mizer model.

Along with these setter functions there are accessor functions for getting the
parameter arrays: `getPredationKernel()`, `getSearchVolume()`, 
`getInteraction()`, `getMaxIntakeRate()`, `getMetabolicRate()`, 
`getExtMortality()`, `getMaturityProportion()`, `getReproductionProportion()`,
`getCatchability()`, `getSelectivity()`, `getPlanktonBirthRate()`,
`getPlanktonCarryingCapacity()`, `getPlanktonParams()`, `getPlanktonDynamics()`,

* Setting of the maximum reproduction rate has been separated out into new
  function `setRmax()`.

## Initial Values and steady state

The MizerParams object now also contains the initial values for the size
spectra. This is particularly useful if the model has been tuned to produce
the observed steady state. The new function `steady()` finds a steady state
for a model and sets it as the initial value. The initial values can be
accessed and changed via functions `initial_n()` and `initial_n_pp()`. The
initial values can be set to the final values of a previous simulation with
`setInitialValues()`.

The MizerParams object now has a slot `initial_effort` that specifies the
  initial fishing effort to which the steady state has been calibrated.

## Extension mechanisms

Mizer now has an extension mechanism that allows other R packages to be
written to generalise the mizer model. More documentation to follow.

## Plotting

* Every plot function now has a plotly version that makes the plot interactive 
  using the plotly package. So for example there is `plotlyBiomass()` as the 
  plotly version of `plotBiomass()`, and so on.
* New `plotGrowthCurves()` plots growth curves and compares them to the von
  Bertalanffy growth curve.
* New `plotDiet()` plots the diet composition as a function of predator size.
* New `highlight` argument to all plot functions that display curves for 
  multiple species. Displays highlighted species with wider lines.
* In the legends of all plots the species are now consistently ordered in the
  same way as in the species parameter data frame.
* All plot functions that are not time-resolved now accept also a MizerParams
  object as an alternative to the MizerSim object to plot the initial state.
* New `plot()` method for MizerParams object to plot the initial state.
* Avoiding duplicate graphs in rmarkdown documents.
* New argument `include_critical` in `plotFeedingLevel()` allows to show also
  the critical feeding level.
* New `wlim` argument to `plotSpectra()` in analogy to the existing `ylim`
  argument to limit the w range in the plot.
* The colours used in plot functions can be set with `setColours()`.
* The default line type is `solid` but this can be changed via the 
  `setLinetypes()` function.
* Use colour and linetype for plots irrespective of the number of species.

## General predation kernel

* Users can now replace the lognormal function in the predation kernel by a
  function of their choice, allowing a differently shaped kernel for each 
  species.
* New `box_pred_kernel()` implements a box-shaped kernel as an alternative to
  the default `lognormal_pred_kernel()`.
* New `power_law_pred_kernel()` implements a power-law kernel with sigmoidal
  cutoffs at both ends. This is suitable for filter feeders.
* Users can sets a predation kernel that has a predator-size-dependent
  predator/prey mass ration (via `setPredationKernel()`). Mizer automatically
  falls back on the old non-FFT code to handle this. (#41)
* New `getPredationKernel()` returns the full 3-dimensional predation kernel array,
  even when this is not stored in MizerParams object.
  
## Other new functions

* New `upgradeParams()` can upgrade MizerParams objects from previous versions 
  of mizer so they work with the new version.
* New `getDiet()` calculates the diet of predators. (#43)
* Alternative functions `RickerRDD()` and `SheperdRDD()` for density-dependence 
  in reproduction, as well as `noRDD()` and `constantRDD()`.
* New gear selectivity function `double_sigmoid_length()` allows modelling
  of escape of large individuals.
* New gear selectivity function `sigmoidal_weight()` is weight-based trawl 
  selectivity function. (Ken H Andersen)
* New `getGrowthCurves()` calculates the growth curves (size at age).
* New `mizerRates()` calculates all the rates needed in the model and collects
  them in a list.
* A convenience function `times()` to extract the times at which simulation 
  results are saved in a MizerSim object.
* Convenience functions `final_n()`, `final_n_pp()` and `final_n_other` to
  access the values at the final time of a simulation.
* New function `getCriticalFeedingLevel()` returns the critical feeding level
  for each species at each size.
* Mizer reexports the `melt()` function from the reshape2 package which allows
  users to convert the arrays returned by mizer functions into data frames
  that can be used for example in ggplot2 and plotly.
* `validSpeciesParams()` checks validity of species parameter data frame and
  sets defaults for missing but required parameters.

## Other new features

* The allometric exponents `n`, `p` and `q` as well as the feeding level `f0`
  can now be set at the species level via columns in `species_params`.
* The critical feeding level `fc` can now be specified as a species parameter 
  and will be used to calculate the metabolic rate parameter `ks` if it is not
  supplied.
* `project()` now shows a progress bar while a simulation is running. Can be
  turned off with `progress_bar = FALSE` argument.
* Satiation can be switched off by setting the maximum intake rate to `Inf`.
* Users can now set their own plankton dynamics instead of the default
  `plankton_semichemostat()`.
* Different species can interact with plankton with different strengths, or not
  feed on plankton at all, controlled by an `interaction_p` column in the
  species parameter data frame.
* The steepness of the maturity ogive can now be controlled via a `w_mat25`
  column in the species parameter dataframe, which gives the size at which
  25% of the individuals of a species are mature.
* The scaling exponent for the allocation of resources into reproduction can
  now be set via the `m` column in the species parameter data frame.
* `project()` can now continue projection from last time step of a previous
  simulation if the first argument is a MizerSim object. The new `append` 
  argument then controls whether the new results are appended to the old.
* Values for minimum plankton size, and minimum and maximum consumer sizes are
  set automatically if not provided in `newMultispeciesParams()`.
* Default values for species parameters are used for missing values within a 
  column in the species parameter data frame, not only if the column is missing 
  entirely.
* Rate functions take defaults for their `n`, `n_pp` and `n_other` arguments
  from the initial values in the `params` argument.
* New `perfect_scaling` argument allows `newTraitParams()` to produce a perfectly 
  scale-invariant model.
* A new `ext_mort_prop` argument in `newTraitParams()` allows the inclusion of
  external mortality.
* Added a data file`NS_params` with the North Sea model MizerParams object.
* Comments can be added to MizerParams objects and any of their slots. Slots
  that have comments are protected from being overwritten with allometric
  defaults.
  
## Documentation

* Mizer now has a documentation website at <https://sizespectrum.org/mizer/>
  for the latest released version and at <https://sizespectrum.org/mizer/dev>
  for the development version. (#48)
* The help pages of mizer functions has been extended massively, see for
  example the help for `newMultispeciesParams()`.
* The vignette chapters are shown as pages on the website.
* The html help pages for plotting functions now show example plots.
* Clarified that mizer uses grams and years as size and time units and is 
  agnostic about whether abundances are per area, per volume or per study area.
  (#42)
* Added a tutorial on using ggplot2 and plotly with mizer.
* Added a tutorial on working with git and GitHub for mizer development.
* Added a FAQ page for developers.
* Added a unit test to automatically run a spell check on documentation.
* Renamed some functions for consistency and to make them easier to understand,
  but kept old names as aliases for backwards compatibility:
  + `getmM2()` -> `getPredMort()`
  + `plotM2` -> `plotPredMort()`
  + `getM2background()` -> `getPlanktonMort()`
  + `getZ()` -> `getMort()`
  + `getESpawning()` -> `getERepro()`
  + `MizerParams()` -> `emptyParams()` or `set_multispecies_model()`
* Renamed maximum reproductive rate from `r_max` to `R_max`.
* Updated list of publications (@Kenhasteandersen)

## Ecosystems
Added ecosystems from N.S. Jacobsen, M. Burgess and K.H. Andersen (2017): Efficiency of fisheries is increasing at the ecosystem level. Fish and Fisheries 18(2) 199- 211. doi:10.1111/faf.12171:
* `data(Benguela_params)` with five species: Anchovy, Sardine, Kingklip, 
  Shallow water hake, and Deep water hake.
* `data(Baltic_params)` with three species: sprat, herring, and cod.
* `data(Barents_params)` with six species.
* `data(NEUSCS_params)` with 24 species.
* `data(NorthSea_params)` with 10 species.

## Bug fixes

* In `getSSB()`, the calculation of the spawning stock biomass is done correctly
  using the maturity ogive instead of the proportion of energy allocated to
  reproduction. (#47)
* The fast FFT method and the old method for calculating integrals now give 
  the same numerical results. (#39)
* `getEncounter()` and `getPredRate()` now set names on the returned arrays.
* Plankton carrying capacity for scale-invariant model is calculated in a way 
  that reduces rounding errors.
* Avoids potential problems with negative numbers due to numerical errors.
* Consistently cutting off predation kernel at 0 and beta + 3 sigma.
* The `ylim` argument is not handled correctly in plots.
* `display_frame()` is now exported.
* plotGrowthCurves() and getGrowthCurves() also works when there is only a 
  single species
* t_start argument in project() is used correctly
* times are not truncated at 3 significant figures, because that would not allow
  something like 2019.
* get_initial_n() gets values for n and q from params object
* `summary()` of MizerParams object reflects the number of non-empty plankton 
  bins. (@patricksykes)
  
## Under the hood

* Increased regression test coverage from 67% to 86%.
* Now using vdiffr package to test plots.
* Converted all S4 methods to functions to decrease the learning curve for
  new developers.
* The calculation of defaults is now handled by new `get_gamma_default()`,
  `get_h_default()` and `get_ks_default()`, making it easier to change or
  extend these in the future.
* Helper function `set_species_param_default()` makes it easier to set default
  values for species parameters.
* Simplified FFT calculations are more readable.
* Using `@inherit` functionality of roxygen2 to reduce duplication in
  roxygen documentation.
* Using `@family' to group function documentation pages.
* The helper functions are now documented and exported.
* `getPhiPrey()` is replaced by `getEncounter()` which now returns the full
  encounter rate, including the contribution from other components. Even
  in the absence of other components, `getEncounter()` differs from the
  old `getPhiPrey()` because it includes the search volume factor.
* Changed naming convention: user-facing function names are now in camelCase.
* Consistently use `params` to refer to an argument of class MizerParams, `sim`
  to refer to an argument of class MizerSim, and `object` to an argument that
  can be either.
* Updated the calls to `setClass()` to follow the new guidelines, replacing
  `representation` by `class` and removing `prototype` and `validity`.
* Added numerical tests.
* Using assert_that to check arguments to functions more often.
* Argument `shiny_progress` renamed to `progress_bar` because they control
  any type of progress bar.
* In documentation renamed "background" to "plankton".
* Using `outer()` instead of `tapply()` where possible to improve readability.
* Avoiding use of `hasArg()` and `anyNA()` because they were not available in R 3.1
* A more robust code for setting up the size grids.
* Improved consistency of when to issue warnings and when to issue messages.
* Changes to MizerParams class:
  + Merged `std_metab` and `activity` slots into a single `metab` slot.
  + Moved `w_min_idx` out of `species_params` into its own slot.
  + Added slot `maturity` to hold the maturity ogive.
  + Added slot `pred_kernel` to hold predation kernel if it has variable
    predator/prey ratio.
  + Added slot `plankton_dynamics` to allow user to specify alternative
    plankton dynamics.
  + Instead of the function in the slot `@srr` we now have the name of the 
    function in `@rate_funcs$RDD`, see #91.
  + Added slots `@other_dynamics`, `@other_params`, `@other_encounter`,
    `@other_mort` and `initial_n_other` to allow mizer extensions to add more 
    ecosystem components.
  + Added slot `@rates_funcs` to allow mizer extensions to replace mizer rate
    functions with their own rate functions.



# mizer 1.0.1

* Now compatible with older versions of R > 3.1.0.
* Skipping a test on CRAN that fails on some machines with different precision.
* Fixing minor typos in documentation.


# mizer 1.0

* Fixed bugs in how the start time of a simulation was handled. This leads to
  small corrections, so that the output of this version is slightly different 
  from previous versions.
* Introduced a scale-invariant trait-based model, set up with 
  `set_scaling_model()`, see section 12 in the vignette.
* Added a function that adds news species to a scale-invariant background, 
  and computes an approximately steady state close to the power law, see
  section 13 in the vignette.
* Created an example shiny app to allow people to use mizer through a web 
  browser without having to install mizer. The app explores the effect of more 
  selective fishing gear in a case study.
* Improvements to plots:
  + Added units to axes
  + Added function for plotting growth curves
  + `PlotYield()` no longer fails when species names are numbers or when a 
     species abundance is zero
  + Added a `total` parameter to several plot functions to add the curve for the 
     total community (sum over all species and plankton)
  + Added a `species` parameter to all plot functions to allow for only a 
      selection of species to be plotted
  + Allow the number of ticks on y-axis in biomass plot to be controlled
* Allow for size- and species-dependent background death.
* Add `@initial_n` and `@initial_n_pp` slots to MizerParams class.
* Now checking that effort times are increasing.
* Corrections in the documentation.
* Improvements to the vignette.
* Add a test of the numeric solution against an analytic solution.


# mizer 0.4

* Improvements made to the speed by evaluating convolution sums via fft,
  removing the bottlenecks in `getPhiPrey()` and `getPredRate()`.
* Using C++ for the inner loop in the project method for extra speed.
* Minor corrections to vignette and documentation to bring them into alignment
  and to document the new home on GitHub and new maintainers.


# mizer 0.3

* Improvements made to the speed of the simulations. Remaining bottle necks 
  are the sweep statements in `getPhiPrey()` and `getPredRate()`.
* Moved tests to new suggested folder.
* Minor changes to documentation to pass new check requirements.


# mizer 0.2

* Release to coincide with the submission of the MEE paper. No major changes. 
  Just minor bug fixes.


# mizer 0.1

* Beta release - just about works but still some gremlins to sort out.
  There are a number of features I'd like to add in the coming releases.
