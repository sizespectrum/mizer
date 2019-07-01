
# mizer 1.0.1.9000 

## Breaking changes

* The default for `min_w_pp`, the smallest size of plankton, has changed in
  `MizerParams()`, `set_trait_model()` and `set_community_model()`.
  It used to be `w_min_pp = 1e-10`. Now the default is the smallest size at
  which any of the species feeds. If your code relies on the old default, you
  now need to supply the `w_min_pp = 1e-10` argument explicitly.
* The default for `max_w`, the largest size in the model, has changed in
  `MizerParams()` and `set_trait_model()`. It used to be 1.1 times the largest 
  asymptotic size of any species in the model otherwise. The unnecessary factor
  of 1.1 has now been removed. If your code relies on the old default, you now
  need to set `w_max` explicitly.
* Removed the `print_it` argument from plot functions.
* plotFeedingLevel() now only plots the values within the size range of each
  species. If for some reason you want the old plots that show a feeding level
  also for sizes that the fish can never have, you need to supply an argument
  `all.sizes = TRUE`.

## Modelling unstructured resources
Besides the size-structured planktonic resource, mizer can now also model any
number of unstructured resource components. Such unstructured components are
appropriate whenever the predation on these components is not size based.
Possible applications include the modelling of detritus as a resource for
detritivores, carrion as a resource for scavengers, or macroflora on which fish
can graze. (#46)

* Each unstructured resource component can have its own dynamics, with its own
  parameters. These are set up with `setResourceDynamics()`.
* Each species can be set to consume resources. The rates are set up with
 `setResourceEncounter()`.
* Example dynamics are provided in `carrion_dynamics()` and `detritus_dynamics()`.

## Setting model parameters
After setting up a mizer model, it is possible to change specific model
parameters with the new functions

* `setPredKernel()`
* `setSearchVolume()`
* `setInteraction()`
* `setIntakeMax()`
* `setMetab()`
* `setBMort()`
* `setReproduction()`
* `setFishing()`
* `setPlankton()`
* `setResourceDynamics()`
* `setResourceEncounter()`

The new function `setParams()` is a wrapper for all of the above functions
and is also used when setting up a new model with `set_multispecies_model()`.
(#51)

## Plotting

* Every plot function now has a plotly version that makes the plot interactive 
  using the plotly package. So for example there is `plotlyBiomass()` as the 
  plotly version of `plotBiomass()`, and so on.
* New `plotGrowthCurves()` plots growth curves and compares them to the von
  Bertalanffy growth curve.
* New `highlight` argument to all plot functions that display curves for 
  multiple species. Displays highlighted species with wider lines.
* In the legends of all plots the species are now consistently ordered in the
  same way as in the species parameter data frame.
* All plot functions that are not time-resolved now accept also a MizerParams
  object as an alternative to the MizerSim object to plot the initial state.
* New `plot()` method for MizerParams object to plot the initial state.
* Avoiding duplicate graphs in rmarkdown documents by setting the default for
  the `print_it` argument in plot functions to `FALSE`.
* New `wlim` argument to `plotSpectra()` in analogy to the existing `ylim`
  argument to limit the w range in the plot.
* New arguments `xlab` and `ylab` for `displayFrames()`.
* Use colour and linetype for plots irrespective of the number of species.
* Plot background species in the colour specified in the `linecolour` slot.

## General feeding kernel

* Users can now replace the lognormal function in the predation kernel by a
  function of their choice, allowing a differently shaped kernel for each 
  species.
* New `box_pred_kernel()` implements a box-shaped kernel as an alternative to
  the default `lognormal_pred_kernel()`.
* Users can sets a predation kernel that has a predator-size-dependent
  predator/prey mass ration (via `setPredKernel()`). Mizer automatically
  falls back on the old non-FFT code to handle this. (#41)
* New `getPredKernel()` returns the full 3-dimensional predation kernel array,
  even when this is not stored in MizerParams object.

## Other new features

* `project()` now shows a progress bar while a simulation is running. Can be
  turned off with `progress_bar = FALSE` argument.
* New `getDiet()` calculates the diet of predators. (#43)
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
* New `getGrowthCurves()` calculates the growth curves (size at age).
* Values for minimum plankton size, and minimum and maximum consumer sizes are
  set automatically if not provided in `set_multispecies_model()`.
* Default values for species parameters are used for missing values within a 
  column in the species parameter data frame, not only if the column is missing 
  entirely.
* New `getRates()` calculates all the rates needed in the model and collects
  them in a list.
* Can set initial state with `setInitial()`.
* Rate functions take defaults for their `initial_n`, `initial_n_pp` and
  `initial_B` arguments from the corresponding slot in the `params` argument.
* Can remove a species from a model with `removeSpecies()`.
* New `perfect` argument allows `set_scaling_model()` to produce a perfectly 
  scale-invariant model.
* In `retuneAbundance()`, Singular Value Decomposition is used to retune 
  background species to keep the community close to power law. Background 
  species with very low abundances are automatically removed.
* Renamed some functions for consistency and to make them easier to understand,
  but kept old names as aliases for backwards compatibility:
  + `getmM2()` -> `getPredMort()`
  + `plotM2` -> `plotPredMort()`
  + `getM2background()` -> `getPlanktonMort()`
  + `getZ()` -> `getMort()`
  + `getESpawning()` -> `getERepro()`
  + `MizerParams()` -> `emptyParams()` or `set_multispecies_model()`
  
## Documentation

* Mizer now has a documentation website at <https://sizespectrum.org/mizer/>
  for the latest released version and at <https://sizespectrum.org/mizer/dev>
  for the development version. (#48)
* The help pages of mizer functions has been extended massively, see for
  example the help for `set_multispecies_model()`.
* The vignette chapters are shown as pages on the website.
* The html help pages for plotting functions now show example plots.
* Clarified that mizer uses grams and years as size and time units and is 
  agnostic about whether abundances are per area, per volume or per study area.
  (#42)

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
  
## Under the hood

* Increased regression test coverage from 67% to 86%.
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
  encounter rate, including the contribution from unstructured resources. Even
  in the absence of unstructured resources, `getEncounter()` differs from the
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
  + Added slot `rho` to hold the resource encounter rates.
  + Added slot `pred_kernel` to hold predation kernel if it has variable
    predator/prey ratio.
  + Added slot `plankton_dynamics` to allow user to specify alternative
    plankton dynamics.
  + Added slots `resource_dynamics` and `resource_params`.
  + Added slot `initial_B` for the initial biomasses of the resources.
* Changes to MizerSim class:
  + Added slot `B` to hold resource biomasses over time


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
* Add initial_n and initial_n_pp slots to mizer params.
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
