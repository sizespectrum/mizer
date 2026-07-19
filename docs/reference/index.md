# Package index

## Overview: the mizer workflow

mizer builds and simulates dynamic, size-structured models of fish
communities. Building and using a model follows five stages, and the
reference sections below are organised in roughly this order:

1.  **Create** a model from species and gear parameters, starting with
    [`newMultispeciesParams()`](https://sizespectrum.org/mizer/reference/newMultispeciesParams.md)
    or one of the simpler
    [`newCommunityParams()`](https://sizespectrum.org/mizer/reference/newCommunityParams.md),
    [`newTraitParams()`](https://sizespectrum.org/mizer/reference/newTraitParams.md)
    and
    [`newSingleSpeciesParams()`](https://sizespectrum.org/mizer/reference/newSingleSpeciesParams.md).
2.  **Calibrate the steady state** so that growth, biomass and yield
    match observations, with
    [`matchGrowth()`](https://sizespectrum.org/mizer/reference/matchGrowth.md),
    [`calibrateBiomass()`](https://sizespectrum.org/mizer/reference/calibrateBiomass.md),
    [`matchBiomasses()`](https://sizespectrum.org/mizer/reference/matchBiomasses.md)
    and
    [`steady()`](https://sizespectrum.org/mizer/reference/steady.md).
3.  **Tune the dynamics** so the model responds realistically to
    perturbations away from the steady state, with
    [`setBevertonHolt()`](https://sizespectrum.org/mizer/reference/setBevertonHolt.md)
    and
    [`setResource()`](https://sizespectrum.org/mizer/reference/setResource.md).
4.  **Project** the model forward in time under a fishing scenario, with
    [`project()`](https://sizespectrum.org/mizer/reference/project.md).
5.  **Analyse and plot** the results, with summary functions such as
    [`getBiomass()`](https://sizespectrum.org/mizer/reference/getBiomass.md)
    and
    [`getYield()`](https://sizespectrum.org/mizer/reference/getYield.md)
    and plots such as
    [`plotBiomass()`](https://sizespectrum.org/mizer/reference/plotBiomass.md)
    and
    [`plotSpectra()`](https://sizespectrum.org/mizer/reference/plotSpectra.md).

New users should start with the [Get started
guide](https://sizespectrum.org/mizer/articles/mizer.html) and the topic
[cheat sheets](https://sizespectrum.org/mizer/articles/index.html), or
open the package overview page below.

- [`mizer`](https://sizespectrum.org/mizer/reference/mizer-package.md)
  [`mizer-package`](https://sizespectrum.org/mizer/reference/mizer-package.md)
  : mizer: Multi-species size-based modelling in R

## Creating a new model

Mizer allows the easy set-up of four different types of models, of
increasing level of complexity. See
<https://sizespectrum.org/mizer/articles/mizer.html#size-spectrum-models>
for a description of these model types. The [Model setup and calibration
cheatsheet](https://sizespectrum.org/mizer/articles/cheatsheet-model-setup-and-calibration.html)
gives a quick overview of the workflow.

- [`newSingleSpeciesParams()`](https://sizespectrum.org/mizer/reference/newSingleSpeciesParams.md)
  **\[experimental\]** : Set up parameters for a single species in a
  power-law background
- [`newCommunityParams()`](https://sizespectrum.org/mizer/reference/newCommunityParams.md)
  : Set up parameters for a community-type model
- [`newTraitParams()`](https://sizespectrum.org/mizer/reference/newTraitParams.md)
  : Set up parameters for a trait-based multispecies model
- [`newMultispeciesParams()`](https://sizespectrum.org/mizer/reference/newMultispeciesParams.md)
  : Set up parameters for a general multispecies model

## Changing model parameters

After you have created a model, you will want to make changes to it
while tuning the model and for investigating the impact of changes in
parameters. See the [Changing model parameters
cheatsheet](https://sizespectrum.org/mizer/articles/cheatsheet-changing-parameters.html)
for a quick reference.

- [`species_params()`](https://sizespectrum.org/mizer/reference/species_params.md)
  [`` `species_params<-`() ``](https://sizespectrum.org/mizer/reference/species_params.md)
  [`is.species_params()`](https://sizespectrum.org/mizer/reference/species_params.md)
  [`given_species_params()`](https://sizespectrum.org/mizer/reference/species_params.md)
  [`is.given_species_params()`](https://sizespectrum.org/mizer/reference/species_params.md)
  [`` `given_species_params<-`() ``](https://sizespectrum.org/mizer/reference/species_params.md)
  [`calculated_species_params()`](https://sizespectrum.org/mizer/reference/species_params.md)
  : Species parameters
- [`gear_params()`](https://sizespectrum.org/mizer/reference/gear_params.md)
  [`` `gear_params<-`() ``](https://sizespectrum.org/mizer/reference/gear_params.md)
  [`is.gear_params()`](https://sizespectrum.org/mizer/reference/gear_params.md)
  : Gear parameters
- [`` `initialN<-`() ``](https://sizespectrum.org/mizer/reference/initialN-set.md)
  [`initialN()`](https://sizespectrum.org/mizer/reference/initialN-set.md)
  : Initial values for fish spectra
- [`` `initialNResource<-`() ``](https://sizespectrum.org/mizer/reference/initialNResource-set.md)
  [`initialNResource()`](https://sizespectrum.org/mizer/reference/initialNResource-set.md)
  : Initial value for resource spectrum
- [`initial_effort()`](https://sizespectrum.org/mizer/reference/initial_effort.md)
  [`` `initial_effort<-`() ``](https://sizespectrum.org/mizer/reference/initial_effort.md)
  : Initial fishing effort
- [`addSpecies()`](https://sizespectrum.org/mizer/reference/addSpecies.md)
  : Add new species
- [`removeSpecies()`](https://sizespectrum.org/mizer/reference/removeSpecies.md)
  : Remove species
- [`renameSpecies()`](https://sizespectrum.org/mizer/reference/renameSpecies.md)
  : Rename species
- [`renameGear()`](https://sizespectrum.org/mizer/reference/renameGear.md)
  : Rename gears
- [`adjustSizeGrid()`](https://sizespectrum.org/mizer/reference/adjustSizeGrid.md)
  : Adjust the size grid
- [`markBackground()`](https://sizespectrum.org/mizer/reference/markBackground.md)
  : Designate species as background species
- [`removeBackgroundSpecies()`](https://sizespectrum.org/mizer/reference/removeBackgroundSpecies.md)
  : Remove all background species
- [`use_predation_diffusion()`](https://sizespectrum.org/mizer/reference/use_predation_diffusion.md)
  [`` `use_predation_diffusion<-`() ``](https://sizespectrum.org/mizer/reference/use_predation_diffusion.md)
  : Get or set the use_predation_diffusion flag
- [`second_order_w()`](https://sizespectrum.org/mizer/reference/second_order_w.md)
  [`` `second_order_w<-`() ``](https://sizespectrum.org/mizer/reference/second_order_w.md)
  **\[experimental\]** : Get or set the second_order_w flags

## Steady state tuning

The first task after creating a multi-species model is to tune the model
parameters so that in its steady state the model reproduces average
observed growth rates, abundances and fisheries yields. The [Model setup
and calibration
cheatsheet](https://sizespectrum.org/mizer/articles/cheatsheet-model-setup-and-calibration.html)
walks through this calibration workflow.

Two families of functions rescale abundances to match observations. The
`calibrate...()` functions apply a single overall scaling factor to the
whole model, whereas the `match...()` functions rescale each species
individually. Within each family, the `...Biomass` variant matches
observed biomasses (a `biomass_observed` column in the species
parameters) while the `...Number` variant matches observed numbers (a
`number_observed` column). Use
[`matchGrowth()`](https://sizespectrum.org/mizer/reference/matchGrowth.md)
to match observed von Bertalanffy growth, and the
`plot...ObservedVsModel()` functions to see how well the current model
reproduces the observations.

- [`steady()`](https://sizespectrum.org/mizer/reference/steady.md) : Set
  initial values to a steady state for the model
- [`steadySingleSpecies()`](https://sizespectrum.org/mizer/reference/steadySingleSpecies.md)
  **\[experimental\]** : Set initial abundances to solution of
  steady-state equation with current rates
- [`matchGrowth()`](https://sizespectrum.org/mizer/reference/matchGrowth.md)
  **\[experimental\]** : Adjust model to produce observed growth
- [`plotBiomassObservedVsModel()`](https://sizespectrum.org/mizer/reference/plotBiomassObservedVsModel.md)
  **\[experimental\]** : Plotting observed vs. model biomass data
- [`calibrateBiomass()`](https://sizespectrum.org/mizer/reference/calibrateBiomass.md)
  **\[experimental\]** : Calibrate the model scale to match total
  observed biomass
- [`calibrateNumber()`](https://sizespectrum.org/mizer/reference/calibrateNumber.md)
  **\[experimental\]** : Calibrate the model scale to match total
  observed number
- [`matchBiomasses()`](https://sizespectrum.org/mizer/reference/matchBiomasses.md)
  **\[experimental\]** : Match biomasses to observations
- [`matchNumbers()`](https://sizespectrum.org/mizer/reference/matchNumbers.md)
  **\[experimental\]** : Match numbers to observations
- [`plotYieldObservedVsModel()`](https://sizespectrum.org/mizer/reference/plotYieldObservedVsModel.md)
  **\[experimental\]** : Plotting observed vs. model yields
- [`scaleModel()`](https://sizespectrum.org/mizer/reference/scaleModel.md)
  **\[experimental\]** : Change scale of the model
- [`scaleRates()`](https://sizespectrum.org/mizer/reference/scaleRates.md)
  **\[experimental\]** : Rescale all rates in a mizer model

## Dynamics tuning

After tuning the steady state, you need to tune the sensitivity of the
dynamics to perturbations away from the steady state. The following
functions allow you to change the model without destroying the steady
state.

- [`setBevertonHolt()`](https://sizespectrum.org/mizer/reference/setBevertonHolt.md)
  : Set Beverton-Holt reproduction without changing the steady state
- [`setResource()`](https://sizespectrum.org/mizer/reference/setResource.md)
  [`resource_rate()`](https://sizespectrum.org/mizer/reference/setResource.md)
  [`` `resource_rate<-`() ``](https://sizespectrum.org/mizer/reference/setResource.md)
  [`resource_capacity()`](https://sizespectrum.org/mizer/reference/setResource.md)
  [`` `resource_capacity<-`() ``](https://sizespectrum.org/mizer/reference/setResource.md)
  [`resource_level()`](https://sizespectrum.org/mizer/reference/setResource.md)
  [`` `resource_level<-`() ``](https://sizespectrum.org/mizer/reference/setResource.md)
  [`resource_dynamics()`](https://sizespectrum.org/mizer/reference/setResource.md)
  [`` `resource_dynamics<-`() ``](https://sizespectrum.org/mizer/reference/setResource.md)
  : Set resource dynamics

## Sharing models

Save a model together with its metadata so it can be archived or shared
with other users.

- [`setMetadata()`](https://sizespectrum.org/mizer/reference/setMetadata.md)
  [`getMetadata()`](https://sizespectrum.org/mizer/reference/setMetadata.md)
  : Set metadata for a model
- [`saveParams()`](https://sizespectrum.org/mizer/reference/saveParams.md)
  [`readParams()`](https://sizespectrum.org/mizer/reference/saveParams.md)
  [`saveSim()`](https://sizespectrum.org/mizer/reference/saveParams.md)
  [`readSim()`](https://sizespectrum.org/mizer/reference/saveParams.md)
  : Save and restore mizer objects

## Running simulations

Project a `MizerParams` object forward in time to produce a `MizerSim`
object containing the full time series of size spectra.

- [`project()`](https://sizespectrum.org/mizer/reference/project.md) :
  Project size spectrum forward in time
- [`projectToSteady()`](https://sizespectrum.org/mizer/reference/projectToSteady.md)
  **\[experimental\]** : Project to steady state

## Accessing results

Extract the raw arrays stored in a `MizerSim` object, such as species
and resource size spectra and fishing effort at each saved time step, or
extract the ecosystem state as a `MizerParams` object.

- [`getParams()`](https://sizespectrum.org/mizer/reference/getParams.md)
  : Extract the model state from a simulation
- [`initialParams()`](https://sizespectrum.org/mizer/reference/initialParams.md)
  : Extract the initial state from a simulation
- [`finalParams()`](https://sizespectrum.org/mizer/reference/finalParams.md)
  : Extract the final state from a simulation
- [`N()`](https://sizespectrum.org/mizer/reference/N.md)
  [`NResource()`](https://sizespectrum.org/mizer/reference/N.md) : Time
  series of size spectra
- [`finalN()`](https://sizespectrum.org/mizer/reference/finalN.md)
  [`finalNResource()`](https://sizespectrum.org/mizer/reference/finalN.md)
  [`idxFinalT()`](https://sizespectrum.org/mizer/reference/finalN.md) :
  Size spectra at end of simulation
- [`getEffort()`](https://sizespectrum.org/mizer/reference/getEffort.md)
  : Fishing effort used in simulation
- [`getTimes()`](https://sizespectrum.org/mizer/reference/getTimes.md) :
  Times for which simulation results are available

## Analysing results

Calculate summary quantities from a `MizerSim` object, such as biomass,
yield, growth, and feeding level, averaged or disaggregated over time,
species, or size. The [Analysis and plotting
cheatsheet](https://sizespectrum.org/mizer/articles/cheatsheet-analysis-and-plotting.html)
gives a quick reference to these functions.

- [`summary_functions`](https://sizespectrum.org/mizer/reference/summary_functions.md)
  : Description of summary functions
- [`getBiomass()`](https://sizespectrum.org/mizer/reference/getBiomass.md)
  : Calculate the total biomass of each species within a size range at
  each time step.
- [`getDiet()`](https://sizespectrum.org/mizer/reference/getDiet.md) :
  Get diet of predator at size, resolved by prey species
- [`getGrowthCurves()`](https://sizespectrum.org/mizer/reference/getGrowthCurves.md)
  : Get growth curves giving weight as a function of age
- [`getN()`](https://sizespectrum.org/mizer/reference/getN.md) :
  Calculate the number of individuals within a size range
- [`getSSB()`](https://sizespectrum.org/mizer/reference/getSSB.md) :
  Calculate the SSB of species
- [`getTrophicLevel()`](https://sizespectrum.org/mizer/reference/getTrophicLevel.md)
  **\[experimental\]** : Get trophic level of individuals at size
- [`getTrophicLevelBySpecies()`](https://sizespectrum.org/mizer/reference/getTrophicLevelBySpecies.md)
  **\[experimental\]** : Get mean trophic level of each species
- [`getYield()`](https://sizespectrum.org/mizer/reference/getYield.md) :
  Calculate the rate at which biomass of each species is fished
- [`getYieldGear()`](https://sizespectrum.org/mizer/reference/getYieldGear.md)
  : Calculate the rate at which biomass of each species is fished by
  each gear
- [`getFeedingLevel()`](https://sizespectrum.org/mizer/reference/getFeedingLevel.md)
  : Get feeding level
- [`getCriticalFeedingLevel()`](https://sizespectrum.org/mizer/reference/getCriticalFeedingLevel.md)
  : Get critical feeding level
- [`w()`](https://sizespectrum.org/mizer/reference/w.md)
  [`w_full()`](https://sizespectrum.org/mizer/reference/w.md)
  [`dw()`](https://sizespectrum.org/mizer/reference/w.md)
  [`dw_full()`](https://sizespectrum.org/mizer/reference/w.md) : Size
  bins

## Calculating rates

Calculate instantaneous ecological rates from a `MizerParams` object,
such as encounter rate, predation mortality, or somatic growth rate.

For readers coming from single-species fisheries assessment, mizer’s
fish mortality rates map onto the standard notation as follows:
predation mortality
[`getPredMort()`](https://sizespectrum.org/mizer/reference/getPredMort.md)
is the multi-species analogue of *M2*, external mortality
[`getExtMort()`](https://sizespectrum.org/mizer/reference/setExtMort.md)
is the residual natural mortality not resolved by the model, fishing
mortality
[`getFMort()`](https://sizespectrum.org/mizer/reference/getFMort.md) is
*F*, and the total mortality
[`getMort()`](https://sizespectrum.org/mizer/reference/getMort.md) is
*Z*, the sum of all of these. The older names
[`getM2()`](https://sizespectrum.org/mizer/reference/getM2.md) and
[`getZ()`](https://sizespectrum.org/mizer/reference/getZ.md) are
retained as deprecated aliases for
[`getPredMort()`](https://sizespectrum.org/mizer/reference/getPredMort.md)
and [`getMort()`](https://sizespectrum.org/mizer/reference/getMort.md).
Note that
[`getResourceMort()`](https://sizespectrum.org/mizer/reference/getResourceMort.md)
is different in kind: it is the predation mortality imposed by fish *on
the background resource* spectrum, not a component of fish mortality
(its deprecated alias is
[`getM2Background()`](https://sizespectrum.org/mizer/reference/getM2Background.md)).

- [`getRates()`](https://sizespectrum.org/mizer/reference/getRates.md) :
  Get all rates
- [`getDiffusion()`](https://sizespectrum.org/mizer/reference/getDiffusion.md)
  : Get diffusion rate from predation
- [`getEGrowth()`](https://sizespectrum.org/mizer/reference/getEGrowth.md)
  : Get energy rate available for growth
- [`getERepro()`](https://sizespectrum.org/mizer/reference/getERepro.md)
  : Get energy rate available for reproduction
- [`getEReproAndGrowth()`](https://sizespectrum.org/mizer/reference/getEReproAndGrowth.md)
  : Get energy rate available for reproduction and growth
- [`getEncounter()`](https://sizespectrum.org/mizer/reference/getEncounter.md)
  : Get encounter rate
- [`getFMort()`](https://sizespectrum.org/mizer/reference/getFMort.md) :
  Get the total fishing mortality rate from all fishing gears by time,
  species and size.
- [`getFMortGear()`](https://sizespectrum.org/mizer/reference/getFMortGear.md)
  : Get the fishing mortality by time, gear, species and size
- [`getFeedingLevel()`](https://sizespectrum.org/mizer/reference/getFeedingLevel.md)
  : Get feeding level
- [`getFlux()`](https://sizespectrum.org/mizer/reference/getFlux.md) :
  Get flux into size bins
- [`getFluxGradient()`](https://sizespectrum.org/mizer/reference/getFluxGradient.md)
  **\[experimental\]** : Get flux gradient
- [`getMort()`](https://sizespectrum.org/mizer/reference/getMort.md) :
  Get total mortality rate
- [`getPredMort()`](https://sizespectrum.org/mizer/reference/getPredMort.md)
  : Get total predation mortality rate
- [`getPredRate()`](https://sizespectrum.org/mizer/reference/getPredRate.md)
  : Get predation rate
- [`getRDD()`](https://sizespectrum.org/mizer/reference/getRDD.md) : Get
  density dependent reproduction rate
- [`getRDI()`](https://sizespectrum.org/mizer/reference/getRDI.md) : Get
  density independent rate of egg production
- [`getResourceMort()`](https://sizespectrum.org/mizer/reference/getResourceMort.md)
  : Get predation mortality rate for resource

## Calculating indicators

Calculate ecological indicators from a `MizerSim` object, such as mean
weight, mean maximum weight, and the Large Fish Index.

- [`indicator_functions`](https://sizespectrum.org/mizer/reference/indicator_functions.md)
  : Description of indicator functions
- [`getCommunitySlope()`](https://sizespectrum.org/mizer/reference/getCommunitySlope.md)
  : Calculate the slope of the community abundance
- [`getMeanMaxWeight()`](https://sizespectrum.org/mizer/reference/getMeanMaxWeight.md)
  : Calculate the mean maximum weight of the community
- [`getMeanWeight()`](https://sizespectrum.org/mizer/reference/getMeanWeight.md)
  : Calculate the mean weight of the community
- [`getProportionOfLargeFish()`](https://sizespectrum.org/mizer/reference/getProportionOfLargeFish.md)
  : Calculate the proportion of large fish

## Plotting results

Visualise size spectra, biomass and yield trajectories, growth curves,
and comparisons of model output with observations. See the [Analysis and
plotting
cheatsheet](https://sizespectrum.org/mizer/articles/cheatsheet-analysis-and-plotting.html)
for a quick reference.

Several plots come in related variants. A plain plot such as
[`plotSpectra()`](https://sizespectrum.org/mizer/reference/plotSpectra.md)
shows a single model or simulation. The `...2` variants
([`plotSpectra2()`](https://sizespectrum.org/mizer/reference/plotSpectra2.md),
[`plotCDF2()`](https://sizespectrum.org/mizer/reference/plotCDF2.md))
overlay **two** objects in one figure so you can compare them, and the
`...Relative` variants
([`plotSpectraRelative()`](https://sizespectrum.org/mizer/reference/plotSpectraRelative.md))
show the ratio between two objects. The `...ObservedVsModel` functions
compare model output against observed data. Most `plot...()` functions
have a matching `get...()` accessor that returns the underlying data
frame if you would rather build the plot yourself.

- [`plotting_functions`](https://sizespectrum.org/mizer/reference/plotting_functions.md)
  : Description of the plotting functions

- [`plot`](https://sizespectrum.org/mizer/reference/plot.md) : Plot
  mizer arrays

- [`plotHover()`](https://sizespectrum.org/mizer/reference/plotHover.md)
  : Create a hover-enabled plotly plot from a mizer object

- [`plot2()`](https://sizespectrum.org/mizer/reference/plot2.md) :
  Compare two mizer arrays in a single plot

- [`plotRelative()`](https://sizespectrum.org/mizer/reference/plotRelative.md)
  : Plot relative difference between two mizer arrays

- [`animate()`](https://sizespectrum.org/mizer/reference/animate.md)
  [`animateSpectra()`](https://sizespectrum.org/mizer/reference/animate.md)
  : Animate size-dependent quantities through time

- [`plotSpectra()`](https://sizespectrum.org/mizer/reference/plotSpectra.md)
  : Plot abundance and biomass spectra

- [`plotSpectra2()`](https://sizespectrum.org/mizer/reference/plotSpectra2.md)
  : Compare abundance and biomass spectra from two objects

- [`plotSpectraRelative()`](https://sizespectrum.org/mizer/reference/plotSpectraRelative.md)
  : Plot relative difference between abundance spectra

- [`plotCDF()`](https://sizespectrum.org/mizer/reference/plotCDF.md) :
  Plot cumulative abundance or biomass distributions

- [`plotCDF2()`](https://sizespectrum.org/mizer/reference/plotCDF2.md) :
  Compare cumulative abundance or biomass distributions from two objects

- [`plotBiomass()`](https://sizespectrum.org/mizer/reference/plotBiomass.md)
  : Plot the biomass of species through time

- [`plotPredMort()`](https://sizespectrum.org/mizer/reference/plotPredMort.md)
  : Plot predation mortality rate of each species against size

- [`plotFeedingLevel()`](https://sizespectrum.org/mizer/reference/plotFeedingLevel.md)
  : Plot the feeding level of species by size

- [`plotYield()`](https://sizespectrum.org/mizer/reference/plotYield.md)
  : Plot the total yield of species through time

- [`plotYieldGear()`](https://sizespectrum.org/mizer/reference/plotYieldGear.md)
  : Plot the total yield of each species by gear through time

- [`plotFMort()`](https://sizespectrum.org/mizer/reference/plotFMort.md)
  : Plot total fishing mortality of each species by size

- [`plot(`*`<MizerParams>`*`)`](https://sizespectrum.org/mizer/reference/plotMizerParams.md)
  :

  Summary plot for `MizerParams` objects

- [`plot(`*`<MizerSim>`*`)`](https://sizespectrum.org/mizer/reference/plotMizerSim.md)
  :

  Summary plot for `MizerSim` objects

- [`plotDiet()`](https://sizespectrum.org/mizer/reference/plotDiet.md)
  **\[experimental\]** : Plot diet, resolved by prey species, as
  function of predator at size.

- [`plotGrowthCurves()`](https://sizespectrum.org/mizer/reference/plotGrowthCurves.md)
  **\[experimental\]** : Plot growth curves

- [`addPlot()`](https://sizespectrum.org/mizer/reference/addPlot.md)
  **\[experimental\]** : Add lines to an existing plot

- [`plotBiomassObservedVsModel()`](https://sizespectrum.org/mizer/reference/plotBiomassObservedVsModel.md)
  **\[experimental\]** : Plotting observed vs. model biomass data

- [`plotYieldObservedVsModel()`](https://sizespectrum.org/mizer/reference/plotYieldObservedVsModel.md)
  **\[experimental\]** : Plotting observed vs. model yields

- [`setColours()`](https://sizespectrum.org/mizer/reference/setColours.md)
  [`getColours()`](https://sizespectrum.org/mizer/reference/setColours.md)
  [`setLinetypes()`](https://sizespectrum.org/mizer/reference/setColours.md)
  [`getLinetypes()`](https://sizespectrum.org/mizer/reference/setColours.md)
  **\[experimental\]** : Set line colours and line types to be used in
  mizer plots

## Setting custom rates

You can override the rates mizer calculates from the species parameters
and gear parameters with your own rate arrays.

- [`setParams()`](https://sizespectrum.org/mizer/reference/setParams.md)
  : Set or change any model parameters
- [`setPredKernel()`](https://sizespectrum.org/mizer/reference/setPredKernel.md)
  [`getPredKernel()`](https://sizespectrum.org/mizer/reference/setPredKernel.md)
  [`pred_kernel()`](https://sizespectrum.org/mizer/reference/setPredKernel.md)
  [`` `pred_kernel<-`() ``](https://sizespectrum.org/mizer/reference/setPredKernel.md)
  : Set predation kernel
- [`setSearchVolume()`](https://sizespectrum.org/mizer/reference/setSearchVolume.md)
  [`getSearchVolume()`](https://sizespectrum.org/mizer/reference/setSearchVolume.md)
  [`search_vol()`](https://sizespectrum.org/mizer/reference/setSearchVolume.md)
  [`` `search_vol<-`() ``](https://sizespectrum.org/mizer/reference/setSearchVolume.md)
  : Set search volume
- [`setInteraction()`](https://sizespectrum.org/mizer/reference/setInteraction.md)
  [`interaction_matrix()`](https://sizespectrum.org/mizer/reference/setInteraction.md)
  [`` `interaction_matrix<-`() ``](https://sizespectrum.org/mizer/reference/setInteraction.md)
  : Set species interaction matrix
- [`setMaxIntakeRate()`](https://sizespectrum.org/mizer/reference/setMaxIntakeRate.md)
  [`getMaxIntakeRate()`](https://sizespectrum.org/mizer/reference/setMaxIntakeRate.md)
  [`intake_max()`](https://sizespectrum.org/mizer/reference/setMaxIntakeRate.md)
  [`` `intake_max<-`() ``](https://sizespectrum.org/mizer/reference/setMaxIntakeRate.md)
  : Set maximum intake rate
- [`setMetabolicRate()`](https://sizespectrum.org/mizer/reference/setMetabolicRate.md)
  [`getMetabolicRate()`](https://sizespectrum.org/mizer/reference/setMetabolicRate.md)
  [`metab()`](https://sizespectrum.org/mizer/reference/setMetabolicRate.md)
  [`` `metab<-`() ``](https://sizespectrum.org/mizer/reference/setMetabolicRate.md)
  : Set metabolic rate
- [`setExtDiffusion()`](https://sizespectrum.org/mizer/reference/setExtDiffusion.md)
  [`ext_diffusion()`](https://sizespectrum.org/mizer/reference/setExtDiffusion.md)
  [`` `ext_diffusion<-`() ``](https://sizespectrum.org/mizer/reference/setExtDiffusion.md)
  : Set external diffusion rate
- [`setExtMort()`](https://sizespectrum.org/mizer/reference/setExtMort.md)
  [`getExtMort()`](https://sizespectrum.org/mizer/reference/setExtMort.md)
  [`ext_mort()`](https://sizespectrum.org/mizer/reference/setExtMort.md)
  [`` `ext_mort<-`() ``](https://sizespectrum.org/mizer/reference/setExtMort.md)
  : Set external mortality rate
- [`setExtEncounter()`](https://sizespectrum.org/mizer/reference/setExtEncounter.md)
  [`getExtEncounter()`](https://sizespectrum.org/mizer/reference/setExtEncounter.md)
  [`ext_encounter()`](https://sizespectrum.org/mizer/reference/setExtEncounter.md)
  [`` `ext_encounter<-`() ``](https://sizespectrum.org/mizer/reference/setExtEncounter.md)
  : Set external encounter rate
- [`setReproduction()`](https://sizespectrum.org/mizer/reference/setReproduction.md)
  [`getMaturityProportion()`](https://sizespectrum.org/mizer/reference/setReproduction.md)
  [`maturity()`](https://sizespectrum.org/mizer/reference/setReproduction.md)
  [`` `maturity<-`() ``](https://sizespectrum.org/mizer/reference/setReproduction.md)
  [`getReproductionProportion()`](https://sizespectrum.org/mizer/reference/setReproduction.md)
  [`repro_prop()`](https://sizespectrum.org/mizer/reference/setReproduction.md)
  [`` `repro_prop<-`() ``](https://sizespectrum.org/mizer/reference/setReproduction.md)
  [`psi()`](https://sizespectrum.org/mizer/reference/setReproduction.md)
  : Set reproduction parameters
- [`setFishing()`](https://sizespectrum.org/mizer/reference/setFishing.md)
  [`getCatchability()`](https://sizespectrum.org/mizer/reference/setFishing.md)
  [`catchability()`](https://sizespectrum.org/mizer/reference/setFishing.md)
  [`` `catchability<-`() ``](https://sizespectrum.org/mizer/reference/setFishing.md)
  [`getSelectivity()`](https://sizespectrum.org/mizer/reference/setFishing.md)
  [`selectivity()`](https://sizespectrum.org/mizer/reference/setFishing.md)
  [`` `selectivity<-`() ``](https://sizespectrum.org/mizer/reference/setFishing.md)
  [`getInitialEffort()`](https://sizespectrum.org/mizer/reference/setFishing.md)
  : Set fishing parameters

## Extending Mizer

See [Extending
mizer](https://sizespectrum.org/mizer/articles/extending-mizer.html),
[Using mizer extension
packages](https://sizespectrum.org/mizer/articles/using-extension-packages.html)
and [Creating a mizer extension
package](https://sizespectrum.org/mizer/articles/creating-extension-packages.html)
for more details.

- [`setRateFunction()`](https://sizespectrum.org/mizer/reference/setRateFunction.md)
  [`getRateFunction()`](https://sizespectrum.org/mizer/reference/setRateFunction.md)
  [`other_params()`](https://sizespectrum.org/mizer/reference/setRateFunction.md)
  [`` `other_params<-`() ``](https://sizespectrum.org/mizer/reference/setRateFunction.md)
  : Set own rate function to replace mizer rate function
- [`setComponent()`](https://sizespectrum.org/mizer/reference/setComponent.md)
  [`removeComponent()`](https://sizespectrum.org/mizer/reference/setComponent.md)
  [`getComponent()`](https://sizespectrum.org/mizer/reference/setComponent.md)
  : Add a dynamical ecosystem component
- [`` `initialNOther<-`() ``](https://sizespectrum.org/mizer/reference/initialNOther-set.md)
  [`initialNOther()`](https://sizespectrum.org/mizer/reference/initialNOther-set.md)
  : Initial values for other ecosystem components
- [`NOther()`](https://sizespectrum.org/mizer/reference/NOther.md)
  [`finalNOther()`](https://sizespectrum.org/mizer/reference/NOther.md)
  : Time series of other components
- [`clearExtensionChain()`](https://sizespectrum.org/mizer/reference/clearExtensionChain.md)
  : Clear the registered extension chain
- [`coerceToExtensionClass()`](https://sizespectrum.org/mizer/reference/coerceToExtensionClass.md)
  : Coerce a mizer object to its registered extension class
- [`getRegisteredExtensions()`](https://sizespectrum.org/mizer/reference/getRegisteredExtensions.md)
  : Get the registered mizer extension chain
- [`recordExtension()`](https://sizespectrum.org/mizer/reference/recordExtension.md)
  : Record an extension and its version stamp on a mizer object
- [`registerExtension()`](https://sizespectrum.org/mizer/reference/registerExtension.md)
  : Register a single mizer extension for this R session
- [`registerExtensions()`](https://sizespectrum.org/mizer/reference/registerExtensions.md)
  : Register mizer extensions for this R session
- [`customFunction()`](https://sizespectrum.org/mizer/reference/customFunction.md)
  **\[experimental\]** : Replace a mizer function with a custom version

## Predation kernels

Functions that determine the size preference of predators for prey,
i.e. the probability of a predator of a given size eating prey of a
given size.

- [`box_pred_kernel()`](https://sizespectrum.org/mizer/reference/box_pred_kernel.md)
  : Box predation kernel
- [`lognormal_pred_kernel()`](https://sizespectrum.org/mizer/reference/lognormal_pred_kernel.md)
  : Lognormal predation kernel
- [`power_law_pred_kernel()`](https://sizespectrum.org/mizer/reference/power_law_pred_kernel.md)
  : Power-law predation kernel
- [`truncated_lognormal_pred_kernel()`](https://sizespectrum.org/mizer/reference/truncated_lognormal_pred_kernel.md)
  : Truncated lognormal predation kernel

## Fishing selectivity functions

Functions that determine the size-selectivity of fishing gears, i.e. the
proportion of fish of a given size that are retained by a gear. The
[Fishing
cheatsheet](https://sizespectrum.org/mizer/articles/cheatsheet-fishing.html)
gives a quick reference to setting up gears, selectivity and effort.

- [`double_sigmoid_length()`](https://sizespectrum.org/mizer/reference/double_sigmoid_length.md)
  : Length based double-sigmoid selectivity function
- [`knife_edge()`](https://sizespectrum.org/mizer/reference/knife_edge.md)
  : Weight based knife-edge selectivity function
- [`sigmoid_length()`](https://sizespectrum.org/mizer/reference/sigmoid_length.md)
  : Length based sigmoid selectivity function
- [`sigmoid_weight()`](https://sizespectrum.org/mizer/reference/sigmoid_weight.md)
  : Weight based sigmoidal selectivity function

## Resource dynamics

Functions governing the time evolution of the background resource
spectrum, together with functions for getting and setting resource
parameters.

- [`resource_constant()`](https://sizespectrum.org/mizer/reference/resource_constant.md)
  : Keep resource abundance constant
- [`resource_logistic()`](https://sizespectrum.org/mizer/reference/resource_logistic.md)
  [`balance_resource_logistic()`](https://sizespectrum.org/mizer/reference/resource_logistic.md)
  : Project resource using logistic model
- [`resource_semichemostat()`](https://sizespectrum.org/mizer/reference/resource_semichemostat.md)
  [`balance_resource_semichemostat()`](https://sizespectrum.org/mizer/reference/resource_semichemostat.md)
  : Project resource using semichemostat model
- [`resource_params()`](https://sizespectrum.org/mizer/reference/resource_params.md)
  [`` `resource_params<-`() ``](https://sizespectrum.org/mizer/reference/resource_params.md)
  : Resource parameters

## Reproduction functions

Functions governing the density-dependent relationship between the
energy invested in reproduction and the actual egg production rate.

- [`BevertonHoltRDD()`](https://sizespectrum.org/mizer/reference/BevertonHoltRDD.md)
  : Beverton Holt function to calculate density-dependent reproduction
  rate
- [`RickerRDD()`](https://sizespectrum.org/mizer/reference/RickerRDD.md)
  **\[experimental\]** : Ricker function to calculate density-dependent
  reproduction rate
- [`SheperdRDD()`](https://sizespectrum.org/mizer/reference/SheperdRDD.md)
  **\[experimental\]** : Sheperd function to calculate density-dependent
  reproduction rate
- [`constantEggRDI()`](https://sizespectrum.org/mizer/reference/constantEggRDI.md)
  **\[experimental\]** : Choose egg production to keep egg density
  constant
- [`constantRDD()`](https://sizespectrum.org/mizer/reference/constantRDD.md)
  **\[experimental\]** : Give constant reproduction rate
- [`noRDD()`](https://sizespectrum.org/mizer/reference/noRDD.md) : Give
  density-independent reproduction rate
- [`getReproductionLevel()`](https://sizespectrum.org/mizer/reference/getReproductionLevel.md)
  : Get reproduction level
- [`getRequiredRDD()`](https://sizespectrum.org/mizer/reference/getRequiredRDD.md)
  : Determine reproduction rate needed for initial egg abundance

## Internal rate functions

These functions are used by
[`project()`](https://sizespectrum.org/mizer/reference/project.md) to
calculate instantaneous rates at each time step. You should use the
`get...()` functions instead of the `project...()` functions.

- [`mizerRates()`](https://sizespectrum.org/mizer/reference/mizerRates.md)
  [`projectRates()`](https://sizespectrum.org/mizer/reference/mizerRates.md)
  : Get all rates needed to project standard mizer model
- [`projectDiffusion()`](https://sizespectrum.org/mizer/reference/mizerDiffusion.md)
  [`mizerDiffusion()`](https://sizespectrum.org/mizer/reference/mizerDiffusion.md)
  : Calculate diffusion rate
- [`projectEGrowth()`](https://sizespectrum.org/mizer/reference/mizerEGrowth.md)
  [`mizerEGrowth()`](https://sizespectrum.org/mizer/reference/mizerEGrowth.md)
  : Get energy rate available for growth needed to project standard
  mizer model
- [`projectERepro()`](https://sizespectrum.org/mizer/reference/mizerERepro.md)
  [`mizerERepro()`](https://sizespectrum.org/mizer/reference/mizerERepro.md)
  : Get energy rate available for reproduction needed to project
  standard mizer model
- [`projectEReproAndGrowth()`](https://sizespectrum.org/mizer/reference/mizerEReproAndGrowth.md)
  [`mizerEReproAndGrowth()`](https://sizespectrum.org/mizer/reference/mizerEReproAndGrowth.md)
  : Get energy rate available for reproduction and growth needed to
  project standard mizer model
- [`projectEncounter()`](https://sizespectrum.org/mizer/reference/mizerEncounter.md)
  [`mizerEncounter()`](https://sizespectrum.org/mizer/reference/mizerEncounter.md)
  : Get encounter rate during projection
- [`projectFMort()`](https://sizespectrum.org/mizer/reference/mizerFMort.md)
  [`mizerFMort()`](https://sizespectrum.org/mizer/reference/mizerFMort.md)
  : Get the total fishing mortality rate from all fishing gears
- [`mizerFMortGear()`](https://sizespectrum.org/mizer/reference/mizerFMortGear.md)
  : Get the fishing mortality needed to project standard mizer model
- [`projectFeedingLevel()`](https://sizespectrum.org/mizer/reference/mizerFeedingLevel.md)
  [`mizerFeedingLevel()`](https://sizespectrum.org/mizer/reference/mizerFeedingLevel.md)
  : Get feeding level needed to project standard mizer model
- [`projectMort()`](https://sizespectrum.org/mizer/reference/mizerMort.md)
  [`mizerMort()`](https://sizespectrum.org/mizer/reference/mizerMort.md)
  : Get total mortality rate needed to project standard mizer model
- [`projectPredMort()`](https://sizespectrum.org/mizer/reference/mizerPredMort.md)
  [`mizerPredMort()`](https://sizespectrum.org/mizer/reference/mizerPredMort.md)
  : Get total predation mortality rate needed to project standard mizer
  model
- [`projectPredRate()`](https://sizespectrum.org/mizer/reference/mizerPredRate.md)
  [`mizerPredRate()`](https://sizespectrum.org/mizer/reference/mizerPredRate.md)
  : Get predation rate needed to project standard mizer model
- [`projectRDI()`](https://sizespectrum.org/mizer/reference/mizerRDI.md)
  [`mizerRDI()`](https://sizespectrum.org/mizer/reference/mizerRDI.md) :
  Get density-independent rate of reproduction needed to project
  standard mizer model
- [`projectResourceMort()`](https://sizespectrum.org/mizer/reference/mizerResourceMort.md)
  [`mizerResourceMort()`](https://sizespectrum.org/mizer/reference/mizerResourceMort.md)
  : Get predation mortality rate for resource needed to project standard
  mizer model
- [`projectRDD()`](https://sizespectrum.org/mizer/reference/projectRDD.md)
  : Get density-dependent reproduction rate during projection

## Internal helper functions

Utility functions used internally by mizer that may also be useful for
users building extensions or working with model objects directly.

- [`age_mat()`](https://sizespectrum.org/mizer/reference/age_mat.md) :
  Calculate age at maturity

- [`age_mat_vB()`](https://sizespectrum.org/mizer/reference/age_mat_vB.md)
  : Calculate age at maturity from von Bertalanffy growth parameters

- [`calc_selectivity()`](https://sizespectrum.org/mizer/reference/calc_selectivity.md)
  : Calculate selectivity from gear parameters

- [`constant_other()`](https://sizespectrum.org/mizer/reference/constant_other.md)
  : Helper function to keep other components constant

- [`default_pred_kernel_params()`](https://sizespectrum.org/mizer/reference/default_pred_kernel_params.md)
  : Set defaults for predation kernel parameters

- [`different()`](https://sizespectrum.org/mizer/reference/different.md)
  : Check whether two objects are different

- [`distanceMaxRelRDI()`](https://sizespectrum.org/mizer/reference/distanceMaxRelRDI.md)
  **\[experimental\]** : Measure distance between current and previous
  state in terms of RDI

- [`distanceSSLogN()`](https://sizespectrum.org/mizer/reference/distanceSSLogN.md)
  **\[experimental\]** : Measure distance between current and previous
  state in terms of fish abundances

- [`emptyParams()`](https://sizespectrum.org/mizer/reference/emptyParams.md)
  : Create empty MizerParams object of the right size

- [`get_f0_default()`](https://sizespectrum.org/mizer/reference/get_f0_default.md)
  : Get default value for f0

- [`get_gamma_default()`](https://sizespectrum.org/mizer/reference/get_gamma_default.md)
  : Get default value for gamma

- [`get_h_default()`](https://sizespectrum.org/mizer/reference/get_h_default.md)
  : Get default value for h

- [`get_initial_n()`](https://sizespectrum.org/mizer/reference/get_initial_n.md)
  : Calculate initial population abundances

- [`get_ks_default()`](https://sizespectrum.org/mizer/reference/get_ks_default.md)
  :

  Get default value for `ks`

- [`get_phi()`](https://sizespectrum.org/mizer/reference/get_phi.md) :
  Get values from feeding kernel function

- [`get_size_range_array()`](https://sizespectrum.org/mizer/reference/get_size_range_array.md)
  : Get size range array

- [`get_steady_state_n()`](https://sizespectrum.org/mizer/reference/get_steady_state_n.md)
  : Calculate steady state abundance

- [`get_time_elements()`](https://sizespectrum.org/mizer/reference/get_time_elements.md)
  : Get array indices for a time range in a MizerSim object

- [`l2w()`](https://sizespectrum.org/mizer/reference/l2w.md)
  [`w2l()`](https://sizespectrum.org/mizer/reference/l2w.md) :
  Length-weight conversion

- [`needs_upgrading()`](https://sizespectrum.org/mizer/reference/needs_upgrading.md)
  : Determine whether a MizerParams or MizerSim object needs to be
  upgraded

- [`project_n()`](https://sizespectrum.org/mizer/reference/project_n.md)
  [`project_n_no_diffusion()`](https://sizespectrum.org/mizer/reference/project_n.md)
  : Project values for first time step of Euler method

- [`project_n_2()`](https://sizespectrum.org/mizer/reference/project_n_2.md)
  : Project values with a predictor-corrector method

- [`project_n_tr_bdf2()`](https://sizespectrum.org/mizer/reference/project_n_tr_bdf2.md)
  : Project values with the TR-BDF2 method

- [`project_simple()`](https://sizespectrum.org/mizer/reference/project_simple.md)
  : Project abundances by a given number of time steps into the future

- [`set_species_param_default()`](https://sizespectrum.org/mizer/reference/set_species_param_default.md)
  : Set a species parameter to a default value

- [`validEffortVector()`](https://sizespectrum.org/mizer/reference/validEffortVector.md)
  : Make a valid effort vector

- [`validGearParams()`](https://sizespectrum.org/mizer/reference/validGearParams.md)
  : Check validity of gear parameters and set defaults

- [`validSpeciesParams()`](https://sizespectrum.org/mizer/reference/validSpeciesParams.md)
  [`validGivenSpeciesParams()`](https://sizespectrum.org/mizer/reference/validSpeciesParams.md)
  : Validate species parameter data frame

- [`valid_gears_arg()`](https://sizespectrum.org/mizer/reference/valid_gears_arg.md)
  : Helper function to assure validity of gears argument

- [`valid_species_arg()`](https://sizespectrum.org/mizer/reference/valid_species_arg.md)
  : Helper function to assure validity of species argument

- [`defaults_edition()`](https://sizespectrum.org/mizer/reference/defaults_edition.md)
  : Default editions

## Classes

The S4 and S3 classes used by mizer, together with functions for
constructing, inspecting, comparing, and validating them.

- [`MizerParams-class`](https://sizespectrum.org/mizer/reference/MizerParams-class.md)
  : A class to hold the parameters for a size based model.

- [`summary(`*`<ArraySpeciesBySize>`*`)`](https://sizespectrum.org/mizer/reference/summary.md)
  [`summary(`*`<ArrayTimeBySpecies>`*`)`](https://sizespectrum.org/mizer/reference/summary.md)
  [`summary(`*`<ArrayTimeBySpeciesBySize>`*`)`](https://sizespectrum.org/mizer/reference/summary.md)
  [`summary(`*`<MizerSim>`*`)`](https://sizespectrum.org/mizer/reference/summary.md)
  [`summary(`*`<MizerParams>`*`)`](https://sizespectrum.org/mizer/reference/summary.md)
  : Summarise mizer objects

- [`str(`*`<ArraySpeciesBySize>`*`)`](https://sizespectrum.org/mizer/reference/str.md)
  [`str(`*`<ArrayTimeBySpecies>`*`)`](https://sizespectrum.org/mizer/reference/str.md)
  [`str(`*`<ArrayTimeBySpeciesBySize>`*`)`](https://sizespectrum.org/mizer/reference/str.md)
  [`str(`*`<MizerSim>`*`)`](https://sizespectrum.org/mizer/reference/str.md)
  [`str(`*`<MizerParams>`*`)`](https://sizespectrum.org/mizer/reference/str.md)
  : Display the structure of mizer objects

- [`compareParams()`](https://sizespectrum.org/mizer/reference/compareParams.md)
  : Compare two MizerParams objects and print out differences

- [`validParams()`](https://sizespectrum.org/mizer/reference/validParams.md)
  : Validate MizerParams object and upgrade if necessary

- [`MizerSim-class`](https://sizespectrum.org/mizer/reference/MizerSim-class.md)
  : A class to hold the results of a simulation

- [`getSimParams()`](https://sizespectrum.org/mizer/reference/getSimParams.md)
  : Extract the projection parameters used to produce a simulation

- [`validSim()`](https://sizespectrum.org/mizer/reference/validSim.md) :
  Validate MizerSim object and upgrade if necessary

- [`MizerSim()`](https://sizespectrum.org/mizer/reference/MizerSim.md) :

  Constructor for the `MizerSim` class

- [`print(`*`<ArraySpeciesBySize>`*`)`](https://sizespectrum.org/mizer/reference/print.md)
  [`print(`*`<ArrayTimeBySpecies>`*`)`](https://sizespectrum.org/mizer/reference/print.md)
  [`print(`*`<ArrayTimeBySpeciesBySize>`*`)`](https://sizespectrum.org/mizer/reference/print.md)
  [`print(`*`<summary.ArraySpeciesBySize>`*`)`](https://sizespectrum.org/mizer/reference/print.md)
  [`print(`*`<summary.ArrayTimeBySpecies>`*`)`](https://sizespectrum.org/mizer/reference/print.md)
  [`print(`*`<summary.ArrayTimeBySpeciesBySize>`*`)`](https://sizespectrum.org/mizer/reference/print.md)
  : Print mizer objects

- [`as.data.frame(`*`<ArraySpeciesBySize>`*`)`](https://sizespectrum.org/mizer/reference/as.data.frame.md)
  [`as.data.frame(`*`<ArrayTimeBySpecies>`*`)`](https://sizespectrum.org/mizer/reference/as.data.frame.md)
  [`as.data.frame(`*`<ArrayTimeBySpeciesBySize>`*`)`](https://sizespectrum.org/mizer/reference/as.data.frame.md)
  : Convert mizer arrays to data frames

- [`ArraySpeciesBySize()`](https://sizespectrum.org/mizer/reference/ArraySpeciesBySize.md)
  [`is.ArraySpeciesBySize()`](https://sizespectrum.org/mizer/reference/ArraySpeciesBySize.md)
  : S3 class for species x size rate arrays

- [`ArrayTimeBySpecies()`](https://sizespectrum.org/mizer/reference/ArrayTimeBySpecies.md)
  [`is.ArrayTimeBySpecies()`](https://sizespectrum.org/mizer/reference/ArrayTimeBySpecies.md)
  : S3 class for time x species arrays

- [`ArrayTimeBySpeciesBySize()`](https://sizespectrum.org/mizer/reference/ArrayTimeBySpeciesBySize.md)
  [`is.ArrayTimeBySpeciesBySize()`](https://sizespectrum.org/mizer/reference/ArrayTimeBySpeciesBySize.md)
  : S3 class for time x species x size arrays

- [`ArrayResourceBySize()`](https://sizespectrum.org/mizer/reference/ArrayResourceBySize.md)
  [`is.ArrayResourceBySize()`](https://sizespectrum.org/mizer/reference/ArrayResourceBySize.md)
  : S3 class for resource size spectra

- [`ArrayTimeByResourceBySize()`](https://sizespectrum.org/mizer/reference/ArrayTimeByResourceBySize.md)
  [`is.ArrayTimeByResourceBySize()`](https://sizespectrum.org/mizer/reference/ArrayTimeByResourceBySize.md)
  : S3 class for time x resource-size arrays

## Example parameter sets

More example parameter sets are available via
<https://sizespectrum.org/mizerExamples>

- [`NS_params`](https://sizespectrum.org/mizer/reference/NS_params.md) :
  Example MizerParams object for the North Sea example
- [`NS_species_params`](https://sizespectrum.org/mizer/reference/NS_species_params.md)
  : Example species parameter set based on the North Sea
- [`NS_species_params_gears`](https://sizespectrum.org/mizer/reference/NS_species_params_gears.md)
  : Example species parameter set based on the North Sea with different
  gears
- [`NS_interaction`](https://sizespectrum.org/mizer/reference/NS_interaction.md)
  : Example interaction matrix for the North Sea example
- [`NS_sim`](https://sizespectrum.org/mizer/reference/NS_sim.md) :
  Example MizerSim object for the North Sea example

## Deprecated

These functions are available for backwards compatibility with earlier
versions of mizer

- [`MizerParams()`](https://sizespectrum.org/mizer/reference/MizerParams.md)
  **\[deprecated\]** :

  Alias for
  [`set_multispecies_model()`](https://sizespectrum.org/mizer/reference/set_multispecies_model.md)

- [`calibrateYield()`](https://sizespectrum.org/mizer/reference/calibrateYield.md)
  **\[deprecated\]** : Calibrate the model scale to match total observed
  yield

- [`completeSpeciesParams()`](https://sizespectrum.org/mizer/reference/completeSpeciesParams.md)
  **\[deprecated\]** :

  Alias for
  [`validSpeciesParams()`](https://sizespectrum.org/mizer/reference/validSpeciesParams.md)

- [`expandSizeGrid()`](https://sizespectrum.org/mizer/reference/expandSizeGrid.md)
  **\[deprecated\]** : Expand the size grid

- [`getESpawning()`](https://sizespectrum.org/mizer/reference/getESpawning.md)
  **\[deprecated\]** :

  Alias for
  [`getERepro()`](https://sizespectrum.org/mizer/reference/getERepro.md)

- [`getM2()`](https://sizespectrum.org/mizer/reference/getM2.md)
  **\[deprecated\]** :

  Alias for
  [`getPredMort()`](https://sizespectrum.org/mizer/reference/getPredMort.md)

- [`getM2Background()`](https://sizespectrum.org/mizer/reference/getM2Background.md)
  **\[deprecated\]** :

  Alias for
  [`getResourceMort()`](https://sizespectrum.org/mizer/reference/getResourceMort.md)

- [`getPhiPrey()`](https://sizespectrum.org/mizer/reference/getPhiPrey.md)
  **\[deprecated\]** : Get available energy

- [`getZ()`](https://sizespectrum.org/mizer/reference/getZ.md)
  **\[deprecated\]** :

  Alias for
  [`getMort()`](https://sizespectrum.org/mizer/reference/getMort.md)

- [`inter`](https://sizespectrum.org/mizer/reference/inter.md)
  **\[deprecated\]** :

  Alias for `NS_interaction`

- [`matchYields()`](https://sizespectrum.org/mizer/reference/matchYields.md)
  **\[deprecated\]** : Match yields to observations

- [`plotM2()`](https://sizespectrum.org/mizer/reference/plotM2.md)
  **\[deprecated\]** :

  Alias for
  [`plotPredMort()`](https://sizespectrum.org/mizer/reference/plotPredMort.md)

- [`setInitialValues()`](https://sizespectrum.org/mizer/reference/setInitialValues.md)
  **\[deprecated\]** : Set initial values to values from a simulation

- [`setRmax()`](https://sizespectrum.org/mizer/reference/setRmax.md)
  **\[deprecated\]** :

  Alias for
  [`setBevertonHolt()`](https://sizespectrum.org/mizer/reference/setBevertonHolt.md)

- [`set_community_model()`](https://sizespectrum.org/mizer/reference/set_community_model.md)
  **\[deprecated\]** : Deprecated function for setting up parameters for
  a community-type model

- [`set_multispecies_model()`](https://sizespectrum.org/mizer/reference/set_multispecies_model.md)
  **\[deprecated\]** : Deprecated obsolete function for setting up
  multispecies parameters

- [`set_trait_model()`](https://sizespectrum.org/mizer/reference/set_trait_model.md)
  **\[deprecated\]** : Deprecated function for setting up parameters for
  a trait-based model
