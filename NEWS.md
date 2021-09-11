# mizer 2.3

## New features

* New plots `plotBiomassObservedVsModel()` and `plotYieldObservedVsModel()`
  contributed by @SamikDatta., together with their plotly counterparts.
* New `calibrateBiomass()`, `calibrateYield()` to set the model scale to agree
  with total observed biomass or total observed yield. Uses the new
  `scaleModel()`.
* New `matchBiomasses()` and `matchYields()` will try to adjust the abundances
  of the species to produce the observed biomasses or yields.
  See blog post at https://bit.ly/2YqXESV .
* There are now accessor and replacement functions for rates. So for example
  instead of `params <- setReproduction(params, maturity = my_maturity)` one
  can simply use `maturity(params) <- my_maturity`. These are documented
  together with the setter functions. #213
* New `setMetadata()` to add information to a MizerParams object describing
  the model, for example a title, a description, the author or list of
  authors, a url and a doi. This will be particularly useful for sharing your
  models with others
* New `saveParams()` for saving a MizerParams object to a file and
  `readParams()` for reading it back in. The resulting files can be shared
  with others who want to run your model.
* A MizerParams object now registers the mizer version under which the model was
  last saved. Should the model not be working as expected in the current version
  of mizer, you can go back to the older version under which presumably it was
  working. This helps with the reproducibility of your research.
* A MizerParams object registers the time when it was created and the time it
  was last modified. See `getMetadata()`. This helps you keep track of 
  different versions of your model.
* `steady()` now has a `preserve` argument with possible values `erepro`,
  `R_max` or `reproduction_level` to specify which quantity to preserve.
  This means that one can continue to use `steady()` also
  once one has started to tune the density dependence in reproduction. #208
* Our website is now using the nice new mizer logo designed by Kira Askaroff
  (www.kiraaskaroff.com)
* There is a new mizer extension package 
  [mizerMR](https://sizespectrum.org/mizerMR/)
  allowing you to include multiple resource spectra in your model.

## Small improvements

* The rownames of `gear_params` are now set to "species, gear", so that one
  can access individual entries with for example
  `gear_params(NS_params)["Cod, Otter", "catchability"]`. #212
* The `z0` argument of `setExtMort()` has been deprecated in favour of
  `ext_mort` in order to avoid confusion with the species parameter `z0`.
* `setColours()` and `setLinetypes()` now issue warnings when invalid values
  are given and ignores NAs.
* The experimental `comment` arguments to the setter functions have been
  removed. #214
* The setter functions have a new `reset` argument which, when set to `TRUE`
  will recalculate the rates from the species_, gear_ and resource_params even
  when custom values had been set. #214
* The `species` argument to various functions, which is checked with 
  `valid_species_arg()`, now does not throw an error even when there is no
  valid species included. Only a warning is issued. That means that for
  example `plotSpectra(NS_params, species = list(), total = TRUE)` is now
  allowed.
* `getComponent()` from the mizer extension mechanism now returns NULL when
  asked for a non-existent component instead of giving an error. This gives
  an easy way to check for the existence of a component.
* The example interaction matrix `inter` for the North Sea model now has the
  alternative name `NS_interaction`, with the old name deprecated.
* Species added with `addSpecies()` are now by default given a reproduction
  level of 1/4 instead of 0, because at the low densities at which they are
  introduced there would otherwise not be enough density dependence to 
  stabilise them.
* The size range arguments `min_w`, `max_w`, `min_l` and `max_l` used in some 
  summary functions and processed by `get_size_range_array()` accept vector
  values setting different limits for different species.
* The resource dynamics function is now also passed the `resource_rate` and the
  `resource_capacity` as arguments, which makes it easier to use them in 
  extension packages.
* Species names are now always coerced to strings, even if the user gives them
  as numbers. #202
* There is a new system for informing the user about how defaults were set by
  `newMultispeciesParams()`, #199
* Many improvements in the documentation.
* Many small improvements to code quality and testing.
* Better social media cards, especially for twitter.
* mizer can be run on binder, https://mybinder.org/v2/gh/sizespectrum/mizer/HEAD?urlpath=rstudio

## Bug fixes

* Changing `linecolour` or `linetype` in the species parameters now actually
  changes the linecolours and linetypes as intended.
* Growth curves calculated with `getGrowthCurves()` and plotted with
  `plotGrowthCurves()` are now correct, and no longer extend above the
  asymptotic size.
* `plotGrowthCurves()` with `species_panel = TRUE` now respects the `species`
  argument to only show growth curves for selected species, it works with
  a MizerParams object as well as a MizerSim object, and it shows the panels
  in the correct order. #201
* Reinstated the example .csv files that were missing from the package because
  the vignettes are no longer included.


# mizer 2.2.1

## New functionality

* The `setBevertonHolt()` function has been expanded with more arguments. It
  allows you to change the density dependence in reproduction without changing
  the steady state of your model.
* The new `getReproductionLevel()` function tells you at what proportion of 
  their maximum reproduction rate the species are operating in your model.
* The package now comes with an example MizerSim object `NS_sim` which holds
  a simulation of the North Sea model.
* New function `plotDataFrame()` allows easier creation of plots.


## Bug fixes

* `setInitialValues()` correctly preserves the gear names on the
  initial effort. Thanks to Axel Rossberg.
* `getFMort()` correctly passes the `t` argument on to any custom fishing
  mortality function you may have written.
* The legends in the plots now only show the species that are actually 
  included in the plot.
  
## Other improvements

* Speed improvement in `mizerPredMort()` suggested by Axel Rossberg.
* `plotSpectra()` now only shows those species in the legend that are
  actually contained in the plot.
* Updated tests of plots to use new version of vdiffr package.
* Some improvements to the examples on the help pages.
* Some functions do more thorough tests of their arguments to give more
  useful error messages.
* `initialNOther()` also works with MizerSim object.
* When `projectToSteady()` is called with `effort`, this effort is now also
  stored in the `initial_effort` slot.
* Improvement to `summary()` which is now using `sprintf()` for better
  formatting and also gives the initial_effort.
* Improved documentation of size grid and bins.
* The arguments to `project_simple()` have been given convenient defaults.
* The tooltips in the plotly plots have been cleaned up a bit.
* Species names are now always coerced to strings, even if the user supplies
  numeric names.
* Update to the "A Multi-Species Model of the North Sea" tutorial to use
  `projectToSteady()`.


# mizer 2.2.0

## New functionality

* New function `newSingleSpeciesParams()` for creating a single species in a 
  power-law background.
* New function `animateSpectra()` creating an animated plot of a simulation.
* New functions `addSpecies()`, `removeSpecies()` and `renameSpecies()`.
* The parameters for an ecosystem component added with `setComponent()` can
  now take any form, they no longer have to be a named list.
* New argument `return_data` in the plot's functions allows to return the 
  data frame used for the ggplot instead of the plot.
  
## Breaking changes

* `steady()` no longer switches off the Beverton-Holt density dependence.
  You can do this manually with `setBevertonHolt()` with `R_factor = Inf`.

## Bug fixes

* `getYield()` now also works with density-dependent fishing mortality.
  Thanks to James Roger for discovering the problem.
* The `gamma` argument now is no longer ignored in `newTraitParams()` but
  correctly overrides the `f0` argument. #188
* `getFMort()` again works correctly when called with a MizerSim object.
* `resource_semichemostat()` no longer fail when at some sizes both the 
  resource growth rate and the resource mortality rate are both zero.
* The default for `no_w` in `newTraitParams()` is now always an integer.
* Problems with different machine precision no longer prompts the error
  "The `w_min_idx` should point to the start of the size bin containing 
  the egg size `w_min`".
* `addSpecies() no longer extends grid due to rounding errors.
* If `valid_species_arg()` is called with `species = NULL` and there are no
  background species then it returns `NULL`.

## Documentation

* New tutorial about single-species sizes-spectrum dynamics.
* Improved documentation of `getDiet()` and `plotDiet()`.
* More info on units added to documentation of summary functions.


# mizer 2.1.0

## New functionality

* New function `projectToSteady()` to run the full dynamics to steady state.
* New functions `distanceSSLogN()` and `distanceMaxRelRDI()` to measure 
  distance between two states.
* New function `compareParams()` to compare two MizerParams objects.
* Added `constantEggRDI()` to allow keeping egg densities fixed.
* When setting custom parameter arrays with the setter functions, it is now
  easy for the user to document that via "comment" arguments. #177
* New function `customFunction()` to allow users to overwrite mizer functions.
* Now if the effort is specified as a named vector giving values only for some 
  gears, the effort for the remaining gears is assumed to be zero.
* Added the possibility to see the output of `plotGrowthCurves` as a panel of 
  species with their respective Von Bertalanffy curves

## Breaking changes

* By default, the functions `plotPredMort()` and `plotFMort` will stop 
  displaying mortality values past the species' asymptotic size. The argument     
  `all.sizes` allows you to continue to show these values.

## Bug fixes

* `getFMort()` now passes time argument correctly. #181
* `validEffortArray()` now sets the dimnames correctly. #173

## Code improvements

* Using `lifecycle` package to indicate status of some functions and arguments as
  'experimental' or 'deprecated'.
* Improved error handling in `setFishing()`. #172
* Made use of vdiffr conditional, as required by §1.1.3.1 of
  'Writing R Extensions'.
* Consistent handling of `species` argument in mizer functions, via the new
  `valid_species_arg()` function. #170
* More tests. Test coverage now at 94.71%
* Improved argument checking in `setInitialValues()`
* Throwing error if `min_w_pp` is larger than `min_w`
* Improved documentation of functions for getting fishing mortality.


# mizer 2.0.4

## Bug fixes

* The value of `t` passed to dynamics functions has been corrected.
* `setReproduction()` now correctly sets the the total proportion psi when the 
  maturity proportion is changed.
  
## Enhancements

* The way times are set in `project()` has been simplified. They are now either
  set by the arguments `t_start`, `t_max` and `t_save` or by the dimension names
  of the `effort` array.
* Renamed `setRmax()` to `setBevertonHolt()` and allow it to work on an
  arbitrary MizerParams object. The old name `setRmax()` is still available as
  alias.
* `mizerFMort()` now can also use the abundances and the rates `e_growth` and 
  `pred_mort`. This is useful for example for implementing balanced harvesting.
* A calculation in the numeric scheme has been simplified.
* `gear_params` is allowed to have zero rows.
* In `validGearParams()` the species name is used as gear name in case 
  `gear_name` is NA.
* `validGearParams()` ensures that all required arguments of the 
  selectivity function are supplied and checks validity of species names.
* `species_params()<-` suppresses warnings.
* When `steady()` fails because RDI is zero it gives a meaningful error message.
* `newCommunityParams()` now protects its zero investment in reproduction with
  a comment.
* The default maturity ogive is truncated at proportions smaller than 1e-8.
* A new helper function `valid_species_arg()` checks validity of species 
  selection arguments.
* `upgradeParams()` can now also upgrade old MizerParams objects that do not 
  have a consistent `initial_effort`.
* A new helper function `validParams()` validates a MizerParams object and 
  automatically upgrades it with `upgradeParams()` if necessary.
* Old MizerParams objects are updated automatically when used in plot functions,
  rate functions, summary functions or in `project()` or `steady()`, #163.
* New function `getRates()` to calculates all rates and collects them in a list.
* `steady()` with `return_sim = TRUE` now creates the MizerSim object the same way 
  as `project()`, namely with the original values in the first time slot.
* Added documentation for `species_params()`, `gear_params()` and
  `resource_params()`.
* Numerous small improvements to documentation.


# mizer 2.0.3

## Bug fixes

* Correct handling of shiny progress bar in `project()`.

## Enhancements

* Consistently passing the time argument to the rate functions. This will
  allow extensions to implement time-dependent rates.
* Passing growth and mortality rate to RDI function.
* Simplified the `getRates()` functions by removing the arguments that passed in
  other rates. Instead the required rates are now always calculated within 
  these functions.
* Improved documentation of rate functions and of how to register your own rate 
  functions.
* In `validGearParams()` handle NAs more gracefully and check that there are
  no duplicates.
* Updated hake-mullet selectivity demonstration shiny app.
* Improved user documentation in several places.


# mizer 2.0.2

## Bug fixes

* Time passed to rate functions is now the actual simulation time, not the time 
  elapsed since start of simulation.
* `upgradeParams()` works also on params objects that were created with a
  development version of mizer.
* When upgrading an older params object, `upgradeParams()` does a better job at 
  guessing the value for `w_pp_cutoff`.
* `getFeedingLevel()`, `getPredMort()`, `setInitialValues()` and `steady()` now
  work also when model has extra components.
* The critical feeding level lines are now mentioned in the legend of
  `plotFeedinglevel()` when called with `include_critical = TRUE`, see #162.
* Avoid annoying warnings from dplyr package when `species_params` is a tibble.

## Name changes

* Renamed the functions `params()`, `effort()` and `times()` to `getParams()`,
  `getEffort()` and `getTimes()` to avoid conflicts.

## Minor enhancements

* Some improvements to documentation.
* More unit tests.
* Uses less memory when time step is very small by not creating array with
  effort values at each time step.
* `getDiet()` also takes into account possible contributions by user-defined
  other components.
* In extension mechanism, now the name of a component is also passed to the
  functions implementing dynamics, encounter and mortality.
* If `project()` is called with `t_max < t_save` then `t_save` is automatically
  reduced so that the result at `t_max` will get saved.
* Start showing progress bar earlier during `project()`.
* New helper function `project_simple()` that projects a given number of
  time steps. This might be useful to extension writers.
* The `...` argument to `project()` is passed on to the dynamics and rate
  functions.
* `steady()` runs faster by using `project_simple()`.
* Documentation on mizer website now has a search bar.


# mizer 2.0.1

## Bug fixes

* `summary()` now also works with non-default feeding kernels. #159
* `validSpeciesParams()` no longer fails when `w_mat25` is not specified. #160
* `setInitialValues()` also works in a model with only a single species. #161
* `resource_params()<-` now works and has unit tests.

## Name changes

Some inconsistencies in the choice of names for parameters was removed by
renaming

* `interaction_p` -> `interaction_resource`
* `r_resource` -> `resource_rate`
* `K_resource` -> `resource_capacity`

## Minor enhancements

* New functions `other_params()<-` and `other_params()` for setting and 
  getting other parameters, for example to be used in user-defined rate
  functions.
* `setInitialValues()` also sets `initial_effort`. #157


# mizer 2.0.0 

This is a major new release with many new features, an internal refactoring of
the code and a new extension mechanism. 

## Backwards compatibility

Nevertheless this version of mizer is almost fully backwards compatible with
version 1.0 with the exception of bug fixes and the following breaking changes:

* The previous version of mizer inconsistently truncated the lognormal predation
  kernel when calculating predation but not when calculating encounter. The new
  version never truncates. That leads to very small differences in simulation
  results.
* Removed the `print_it` argument from plot functions.
* `plotFeedingLevel()` now only plots the values within the size range of each
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
  (https://sizespectrum.org/mizerExperimental/)
* During runs of `project()` a progress bar is displayed by default. You can 
  turn this off with the option `progress_bar = FALSE.
* Throughout mizer the term "plankton" has been replaced by "resource", which
  affects the labelling of the resource spectrum in plots.

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

* `species_params<-()`
* `resource_params<-()`
* `gear_params<-()`
* `setPredKernel()`
* `setSearchVolume()`
* `setInteraction()`
* `setMaxIntakeRate()`
* `setMetabolicRate()`
* `setExtMort()`
* `setReproduction()`
* `setFishing()`
* `setResource()`

The new function `setParams()` is a wrapper for all of the above functions
and is also used when setting up a new model with `newMultispeciesParams()`.
(#51)

The documentation for these functions serves to explain the details of the
mizer model.

Along with these setter functions there are accessor functions for getting the
parameter arrays: `getPredKernel()`, `getSearchVolume()`, 
`getInteraction()`, `getMaxIntakeRate()`, `getMetabolicRate()`, 
`getExtMort()`, `getMaturityProportion()`, `getReproductionProportion()`,
`getCatchability()`, `getSelectivity()`, `getResourceRate()`,
`getResourceCapacity()`, `getResourceParams()`, `getResourceDynamics()`,

* Setting of the maximum reproduction rate has been separated out into new
  function `setRmax()`.

## Initial Values and steady state

The MizerParams object now also contains the initial values for the size
spectra. This is particularly useful if the model has been tuned to produce
the observed steady state. The new function `steady()` finds a steady state
for a model and sets it as the initial value. The initial values can be
accessed and changed via functions `initialN()` and `initialNResource()`. The
initial values can be set to the final values of a previous simulation with
`setInitialValues()`.

The MizerParams object now has a slot `initial_effort` that specifies the
  initial fishing effort to which the steady state has been calibrated.

## Extension mechanisms

Mizer now has an extension mechanism that allows other R packages to be
written to generalise the mizer model. See `setRateFunction()` and
`setComponent()`. This mechanism is still experimental and may change as we
gain experience in writing extensions for mizer.

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
* Avoiding duplicate graphs in R Markdown documents.
* New argument `include_critical` in `plotFeedingLevel()` allows to show also
  the critical feeding level.
* New `wlim` argument to `plotSpectra()` in analogy to the existing `ylim`
  argument to limit the `w` range in the plot.
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
  predator/prey mass ration (via `setPredKernel()`). Mizer automatically
  falls back on the old non-FFT code to handle this. (#41)
* New `getPredKernel()` returns the full 3-dimensional predation kernel array,
  even when this is not stored in MizerParams object.
  
## New gear setup
Now it is finally possible to have several gears (or fleets) targeting the same
species. The information is set up via a new `gear_params()` data frame. See
`setFishing()` for details.
  
## Other new functions

* There are now accessor functions for all slots in the MizerParams and
  MizerSim objects. For example to get at the size grid and its spacing you 
  would now use `w()`, `w_full()`, `dw()`, `dw_full()`.
* New `upgradeParams()` and `upgradeSim()` can upgrade objects from 
  previous versions of mizer so they work with the new version.
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
* Convenience functions `finalN()`, `finalNResource()` and `finalNOther()` as
  well as `idxFinalT()` to access the values at the final time of a simulation.
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
* Users can now set their own resource dynamics instead of the default
  `resource_semichemostat()`.
* Different species can interact with resource with different strengths, or not
  feed on resource at all, controlled by an `interaction_resource` column in the
  species parameter data frame.
* The steepness of the maturity ogive can now be controlled via a `w_mat25`
  column in the species parameter dataframe, which gives the size at which
  25% of the individuals of a species are mature.
* The scaling exponent for the allocation of energy into reproduction can
  now be set via the `m` column in the species parameter data frame.
* `project()` can now continue projection from last time step of a previous
  simulation if the first argument is a MizerSim object. The new `append` 
  argument then controls whether the new results are appended to the old.
* Values for minimum resource size, and minimum and maximum consumer sizes are
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
* Gear selectivity functions now can use the species parameters.
  
## Documentation

* Mizer now has a documentation website at <https://sizespectrum.org/mizer/>
  for the latest released version and at <https://sizespectrum.org/mizer/dev/>
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
  + `getM2background()` -> `getResourceMort()`
  + `getZ()` -> `getMort()`
  + `getESpawning()` -> `getERepro()`
  + `MizerParams()` -> `emptyParams()` or `set_multispecies_model()`
* Renamed maximum reproductive rate from `r_max` to `R_max`.
* Updated list of publications (@Kenhasteandersen)
* Using R Markdown in all roxygen comments.

## Bug fixes

* In `getSSB()`, the calculation of the spawning stock biomass is done correctly
  using the maturity ogive instead of the proportion of energy allocated to
  reproduction. (#47)
* The fast FFT method and the old method for calculating integrals now give 
  the same numerical results. (#39)
* `getEncounter()` and `getPredRate()` now set names on the returned arrays.
* Resource carrying capacity for scale-invariant model is calculated in a way 
  that reduces rounding errors.
* Avoids potential problems with negative numbers due to numerical errors.
* Consistently cutting off predation kernel at 0 and beta + 3 sigma.
* The `ylim` argument is not handled correctly in plots.
* `display_frame()` is now exported.
* `plotGrowthCurves()` and `getGrowthCurves()` also works when there is only a 
  single species
* `t_start` argument in `project()` is used correctly
* times are not truncated at 3 significant figures, because that would not allow
  something like 2019.
* `get_initial_n()` gets values for `n` and `q` from params object
* `summary()` of MizerParams object reflects the number of non-empty resource 
  bins. (@patricksykes)
  
## Under the hood

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
* In documentation renamed "background" and "plankton" consistently to "resource".
* Using `outer()` instead of `tapply()` where possible to improve readability.
* Avoiding use of `hasArg()` and `anyNA()` because they were not available in R 3.1
* A more robust code for setting up the size grids.
* Improved consistency of when to issue warnings and when to issue messages.
* Split large code files into smaller files.
* Changes to MizerParams class:
  + Merged `@std_metab` and `@activity` slots into a single `@metab` slot.
  + Moved `@w_min_idx` out of `@species_params` into its own slot.
  + Added slot `@maturity` to hold the maturity ogive.
  + Added slot `@pred_kernel` to hold predation kernel if it has variable
    predator/prey ratio.
  + Added slot `@resource_dynamics` to allow user to specify alternative
    resource dynamics.
  + Added slot `@gear_dynamics` to species to be targeted by multiple gears.
  + Added slot `@ft_mask` that is used when calculating predation rates using
    the Fourier transform method.
  + Added slot `@rates_funcs` to allow mizer extensions to replace mizer rate
    functions with their own rate functions.
  + Instead of the function in the slot `@srr` we now have the name of the 
    function in `@rate_funcs$RDD`, see #91.
  + Added slots `@other_dynamics`, `@other_params`, `@other_encounter`,
    `@other_mort` and `@initial_n_other` to allow mizer extensions to add more 
    ecosystem components.



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
     total community (sum over all species and resource)
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
