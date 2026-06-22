# mizer (development version)

- Extension packages can now upgrade their own data in saved model objects
  independently of the mizer version. The `@extensions` slot can record, for
  each extension, the version of the extension package that the object conforms
  to (write entries with the new `recordExtension()`). `needs_upgrading()` flags
  an object when an extension's recorded version is missing or older than the
  installed package version, and `validParams()` then calls the extension's own
  `upgrade()` method (an S3 method on `utils::upgrade()`, registered with
  `@exportS3Method utils::upgrade`). The core mizer upgrade is now itself the
  `upgrade.MizerParams()` / `upgrade.MizerSim()` method. See the "Upgrading
  objects across versions of your extension" section of
  `vignette("creating-extension-packages")`.
  
- `getTrophicLevel()` and `getTrophicLevelBySpecies()` now assign the resource a
  size-dependent trophic level
  `T_R(w) = max(1, 1 + log(w / w_R) / log(beta_R))` instead of treating it as
  trophic level 0. The new `w_R` (average primary-producer size) and `beta_R`
  (average resource predator/prey mass ratio) arguments control this.

- Resource functions now return classed objects that support the same
  convenient `print()`, `summary()`, `plot()`, and `as.data.frame()` methods
  as the consumer rate functions. `getResourceMort()`, `initialNResource()`,
  `finalNResource()`, `resource_rate()`, `resource_capacity()`, and
  `resource_level()` return an `ArrayResourceBySize` object, and `NResource()`
  returns an `ArrayTimeByResourceBySize` object. So you can now do e.g.
  `plot(getResourceMort(NS_params))` or `plot(NResource(NS_sim))`.
  
- Fixed a bug in `project()` where the abundances of other components (set via
  `setComponent()`) were advanced only once per saved time step instead of once
  per `dt` time step. Their dynamics are now integrated with the same time step
  as the consumer and resource spectra, so results no longer depend on `t_save`.
  Under the second-order time-stepping methods (`"predictor_corrector"` and
  `"tr_bdf2"`) the other components now also receive a corrector step with the
  midpoint rates, so they are integrated to the same second-order accuracy as
  the resource and consumer spectra instead of being left at first order.

- `MizerParams` gains an `ft_pred_kernel_d` slot holding a third Fourier-space
  predation kernel, used by the predation-diffusion integral (when
  `use_predation_diffusion` is `TRUE`). When the `bin_average` entry of
  `second_order_w` is `TRUE`, this kernel carries the extra power of prey size
  that the diffusion integrand needs (`w_p^2 dw_p`, the `β^{3s}` Jacobian)
  instead of reusing the encounter kernel's `β^{2s}`, so the predation-diffusion
  rate is now also second order. In the default first-order scheme it equals
  `ft_pred_kernel_e`, so existing models are byte-identical. (#384)
  Additionally, when `second_order_w[["bin_average"]]` is `TRUE`, this kernel is
  now correctly predator-bin averaged (via a trapezoid fold over adjacent
  offsets), matching the second-order requirements of the diffusion transport step.


- `plotYield()` now uses `sim2 = NULL` instead of `missing(sim2)` to detect
  the optional second simulation argument. This is backward-compatible and
  makes the function work correctly with `do.call()`.
  
- `getRDI()`, `getRDD()`, and `getFlux()` on a `MizerSim` object now correctly
  use the simulated time-varying effort instead of the initial effort. (#370)
  
- Added a new vignette explaining the calculation of default parameter values
  (#189).

- `summary()` of a `MizerSim` object now reports the fishing effort that was
  actually used during the simulation rather than the model's `initial_effort`.
  Gears whose effort varied over time show the mean effort, flagged with a note
  giving the min-max range.
  
- The `second_order_w` slot is now a named list with a character entry `flux`
  (`"upwind"`, `"van_leer"`, or `"centred"`) and a logical entry `bin_average`.
  Setting `second_order_w(params) <- TRUE` enables both second-order features:
  `flux = "van_leer"` (TVD advective reconstruction, keeps abundances
  non-negative) and `bin_average = TRUE` (bin-averaged sinks and quadratures).
  The `"centred"` flux (unlimited, genuinely second order at extrema) can be
  selected with `second_order_w(params) <- c(flux = "centred")`. Old objects
  with the previous logical `second_order_w` slot are upgraded automatically.
  The default first-order upwind scheme is unchanged. See [second_order_w()]
  and the "Numerical Details" vignette.

- When the `bin_average` entry of the `second_order_w` slot is `TRUE`, the
  predation-rate kernel `ft_pred_kernel_p` is averaged over the
  **prey** bin (a trapezoid fold of the kernel), completing the predator-bin
  integral above. Predation mortality is a sink integrated against the prey
  density over the prey bin, so the prey-bin average is the form it needs to be
  second order; `getPredMort()` and `getResourceMort()` (and `getPredRate()`)
  pick this up automatically with no change to the rate functions and no extra
  runtime cost. The default (point-sampled) behaviour is unchanged.

- When the `bin_average` entry of the `second_order_w` slot is `TRUE`,
  `setExtMort()` now replaces the point-sampled power-law external mortality
  \eqn{z_{ext} w^d} by its exact bin average over each bin, making the external
  mortality sink consistent with the finite-volume scheme. The default
  (point sampling) is unchanged.

- When the `bin_average` entry of the `second_order_w` slot is `TRUE`,
  `setExtDiffusion()` likewise replaces the point-sampled power-law external
  diffusion \eqn{D_{ext} w^{n+1}} by its exact bin average, exactly as
  `setExtMort()` does for mortality. The diffusion coefficient is a rate
  represented as a bin average, so whether it is point-sampled or bin-averaged
  is now governed by `bin_average` (the predation-diffusion contribution already
  follows the bin-averaged predation kernel); `flux_limiter` only governs the
  advective reconstruction. A fully second-order transport step therefore needs
  both flags, and the second-order transport routine consumes the diffusion from
  `getDiffusion()` directly rather than averaging it itself. The default
  (point sampling) is unchanged.

- When the `bin_average` entry of the `second_order_w` slot is `TRUE`, the
  reproduction integral in `mizerRDI()` now trapezoidally bin-averages the
  per-capita reproductive investment \eqn{\psi(w) E_r(w)} against the
  cell-average abundance instead of taking its left-bin-edge value, making the
  density-independent reproduction rate second order in the bin size. The full
  per-capita investment is averaged together (not `psi` alone) to capture the
  variation of both factors across the bin. The default (left-edge) behaviour
  is unchanged.
  
- When the `bin_average` entry of the `second_order_w` slot is `TRUE`,
  `setResource()` now builds the auto-calculated resource rate
  \eqn{r_p w^{n-1}} and carrying capacity \eqn{\kappa w^{-\lambda}} from their
  exact bin averages over each size bin instead of point-sampling at the left
  bin edge, with the bin straddling `w_pp_cutoff` getting the partial average.
  This makes the unpredated semichemostat equilibrium equal to the cell-average
  of the background spectrum that the bin-integrated encounter convolution (#374)
  consumes. Only auto-calculated (scalar) rates and capacities are affected;
  user-supplied full vectors are left untouched, and the default
  (point-sampled) behaviour is unchanged. Note that enabling this shifts the
  resource spectrum by `O(Δw)` (more on coarse grids), so calibrated models may
  need recalibrating.

- `project()` gains a new time-stepping option `method = "tr_bdf2"`. This is an
  L-stable, second-order TR-BDF2 scheme that retains the second-order accuracy
  of `method = "predictor_corrector"` while damping the oscillations that the
  Crank-Nicolson corrector can show at large time steps. Like the other methods
  it only requires tridiagonal solves. See the "Numerical Details" vignette.
  
- Under the second-order methods (`"predictor_corrector"` and `"tr_bdf2"`) the
  resource is now advanced with midpoint resource mortality rather than the
  start-of-step value, so the resource spectrum is also second order in time.
  The `"euler"` method and the steady states are unchanged.

- `MizerParams` gains a `second_order_w` slot — a named logical vector with
  entries `flux_limiter` and `bin_average` (both default `FALSE`). When `TRUE`,
  these enable a second-order advective flux and second-order bin-averaged rate
  quadratures respectively. Both need to be `TRUE` for a fully second-order
  scheme. Use the new `second_order_w()` / `second_order_w<-()` accessors to
  get and set the flags. The setter accepts a single logical (sets both entries)
  or a named vector for individual control. It re-runs `setParams()` to rebuild
  precomputed arrays.

- When the `bin_average` entry of `second_order_w` is `TRUE`, the summary
  integrals (`getBiomass()`, `getSSB()`, `getYield()`, `getYieldGear()`,
  `getDiet()`, `getTrophicLevel()`) now use the trapezoidal bin-average of the
  size weight rather than its left-bin-edge value, making these diagnostics
  second order in the bin size. `getN()` is unchanged (its weight is already
  exact). With the default `bin_average = FALSE` the outputs are unchanged.
  Note that enabling bin-averaging shifts reported biomass/yield/SSB by
  `O(Δw)`, so calibrated models may need recalibrating.
  
- When `second_order_w()[["bin_average"]]` is `TRUE`, `calc_selectivity()` now
  stores the bin average of the selectivity over each size bin instead of its
  value at the left bin edge, making the fishing mortality second order in the
  bin size. A knife-edge gear then gets the exact fraction of the straddling bin
  that lies above the knife edge. The default (point-sampled) behaviour is
  unchanged.


# mizer 3.0.0

This release brings new biological realism, improved numerics, a richer 
interactive analysis experience, and a composable extension framework.
For an overview see the
[blog post](https://blog.mizer.sizespectrum.org/posts/2026-05-13-mizer-3-0-announcement/)
pre-announcing the release.

## Diffusion in mizer

The McKendrick-von Foerster equation now supports a diffusion term, allowing
individual variability in growth to be modelled.

- New `getDiffusion()` calculates the total diffusion rate D(w) (g²/year) for
  each species, combining the predation-induced diffusion from the jump-growth
  equation and any externally specified diffusion set via `setExtDiffusion()`.
  It has both `MizerParams` and `MizerSim` methods and returns an
  `ArraySpeciesBySize` or `ArrayTimeBySpeciesBySize` object respectively,
  consistent with the other rate-getter functions.

- The external diffusion coefficient is held in a new `ext_diffusion` slot in
  `MizerParams`. Use `setExtDiffusion()` / `ext_diffusion()` /
  `ext_diffusion<-()` to set and retrieve it. The new species parameter `D_ext`
  (default 0) sets the coefficient of an external diffusion power law;
  `setExtDiffusion()` calculates the default array from species parameters when
  no custom array is supplied, following the same pattern as
  `setExtEncounter()`.

- `MizerParams` gains a `use_predation_diffusion` slot (logical, default
  `FALSE`). When `FALSE` (the default), `mizerDiffusion()` omits the
  predation-induced diffusion term, preserving the behaviour of previous mizer
  versions. Set to `TRUE` via the new `use_predation_diffusion()` accessor to
  enable the jump-growth diffusion term.

- New `getFlux()` function calculates the flux of individuals entering each
  size class, combining the advective flux from somatic growth and the
  diffusive flux. It has a `power` argument, similar to that of `plotSpectra()`,
  for multiplying the flux by a power of the weight; `power = 1` gives the flux
  of biomass.

- `getRequiredRDD()` is exported. It calculates the recruitment rate needed
  to maintain a given initial abundance, accounting for both growth and
  diffusion.

- `steadySingleSpecies()` correctly preserves the steady state under
  `project()`, including when diffusion is non-zero.

- The vignette [cohort dynamics](https://sizespectrum.org/mizer/articles/cohort_dynamics_and_diffusion.html)
  demonstrates the effect of diffusion in an example.


## Higher-order numerical scheme

- `project()`, `projectToSteady()` and `steady()` gain a `method` argument for
  choosing the consumer density time-stepper. The default `"euler"` preserves
  the existing semi-implicit update, while `"predictor_corrector"` uses a new
  second-order predictor-corrector method. The accuracy of the two methods is
  compared in the [numerical details](https://sizespectrum.org/mizer/articles/numerical_details.html)
  vignette.

- `MizerSim` objects now have a `sim_params` slot (a named list) that records
  the projection parameters — currently `method` and `dt` — passed to
  `project()` or `projectToSteady()`. The new `getSimParams()` accessor
  retrieves this list. When `project()` is called on an existing `MizerSim`
  object it defaults `dt` and `method` from the stored `sim_params`, with a
  warning if the supplied values differ. Older objects are upgraded
  automatically by `validSim()`, with `sim_params` set to an empty list.

- `project_n()` and `project_n_2(2)` are new exported functions, factored out of
  `project_simple()`, that projects the abundance spectrum forward in time with
  the different methods. 

## Convenient plot methods for mizer return values

- New `ArraySpeciesBySize` S3 class for the species × size arrays returned by
  many mizer functions. An `ArraySpeciesBySize` object behaves like a regular
  matrix for arithmetic and subsetting but carries a human-readable
  `value_name` and `units` attribute and provides enhanced `print()`,
  `summary()`, `plot()`, and `as.data.frame()` methods. The `plot()` method
  accepts `log_y`, `wlim`, and `ylim` arguments for controlling the y-axis
  scale and limits.

- New `ArrayTimeBySpecies` S3 class for the time × species arrays returned by
  `getBiomass()`, `getSSB()`, `getN()`, and `getYield()` when called on a
  `MizerSim` object. Like `ArraySpeciesBySize`, it carries `value_name` and
  `units` attributes and provides enhanced `print()`, `summary()`, `plot()`,
  and `as.data.frame()` methods. The `plot()` method accepts `log` and `ylim`
  arguments.

- New `ArrayTimeBySpeciesBySize` S3 class for the time × species × size arrays.
  The `N()` accessor on a `MizerSim` object now returns an
  `ArrayTimeBySpeciesBySize` object. Many rate-getter functions —
  `getEGrowth()`, `getEReproAndGrowth()`, `getPredMort()`, `getFMort()`,
  `getMort()`, `getFeedingLevel()`, `getEncounter()`, `getPredRate()`,
  `getRDI()`, `getRDD()` — now also accept a `MizerSim` object and return an
  `ArrayTimeBySpeciesBySize`. An `animate()` method allows interactive
  playback. Subsetting an `ArrayTimeBySpeciesBySize` object returns an
  `ArraySpeciesBySize` object when a single time is selected, and an
  `ArrayTimeBySpecies` object when a single size is selected.

- New `plot2()` generic with methods for comparing two compatible mizer array
  objects in one plot, with species or group shown by colour and model by
  linetype. The `plotSpectra2()` helper has moved from `mizerExperimental` into
  mizer for comparing two abundance spectra.

- New `plotRelative()` generic with methods for plotting the symmetric relative
  difference between two compatible mizer array objects. The
  `plotSpectraRelative()` and `plotlySpectraRelative()` helpers have moved from
  `mizerExperimental` into mizer.

- New `plotCDF()` and `plotCDF2()` generics for plotting cumulative abundance
  or biomass distributions from `MizerParams` and `MizerSim` objects, together
  with `plotlyCDF()` and `plotlyCDF2()` wrappers.

- New `plotHover()` generic with methods for `ArraySpeciesBySize`,
  `ArrayTimeBySpecies`, `ArrayTimeBySpeciesBySize`, and `mizer_plot` converts
  mizer plots into hover-enabled plotly figures.

- New `addPlot()` generic with methods for adding `ArraySpeciesBySize` and
  `ArrayTimeBySpecies` values as extra lines on an existing compatible ggplot.

- The `animate()` methods produces animated plots showing the time evolution
  during a simulation. It can take a`MizerSim` and `ArrayTimeBySpeciesBySize`
  argument and supports axis range settings (`xlim`, `ylim`), timing controls,
  interpolation options, arguments `log_x` `log_y` and `log` to control which
  axis is log-transformed, and `total` and `background` arguments, consistent
  with `plotSpectra()`.

- Plotting functions now consistently expose `log_x`, `log_y` and `log`
  arguments. In all cases, when supplied, `log` overrides `log_x` and `log_y`.
  `plotBiomass()` and `plotYield()` keep support for logical `log` values for
  backward compatibility.

- Time-filtering is now consistent across all time-series plot functions via a
  new `tlim` parameter (analogous to `wlim` and `ylim`): a length-two numeric
  vector `c(start, end)` that restricts the plotted time window. `plotYield()`,
  `plotYieldGear()`, and `animate()` gain this parameter for the first time.
  `plotBiomass()` and `animate.MizerSim()` now use `tlim` in place of the
  former `start_time`/`end_time` and `time_range` parameters respectively;
  the old parameters are deprecated and will be removed in a future release.

- Size-based plots now accept `size_axis = "l"` to show length in cm on the
  size axis instead of weight in grams, using the species' allometric
  weight-length relationship.

- Size-based plots with a `size_axis` argument now accept `llim`, the
  length-axis equivalent of `wlim`, for filtering and limiting plots when
  `size_axis = "l"`.


## Extracting model state from a simulation

- A shift in interpretation of a MizerParams object from just a specification
  of the model to a representation of its state, consisting of both model
  parameters and current values of the state variables (the abundances).

- `getParams(sim, time_range, geometric_mean = FALSE)` now extracts the
  ecosystem state from a `MizerSim` object at a particular time or averaged
  over a time range. When no `time_range` is given, the state at the final time
  step is extracted. New `finalParams(sim)` and `initialParams(sim)` return the
  states at the initial and final times of a simulation respectively.

- Once a state has been extracted from a simulation, it can be analysed by all
  the existing mizer functions. For that purpose the indicator functions
  `getProportionOfLargeFish()`, `getMeanWeight()`, `getMeanMaxWeight()`, and
  `getCommunitySlope()` now also accept a `MizerParams` object and return a
  single value (or named vector for `getMeanMaxWeight()` with 
  `measure = "both"`) calculated from that state. Closes #262.

- `setInitialValues()` is deprecated. Replace
  `setInitialValues(params, sim)` with `finalParams(sim)` (or
  `getParams(sim, time_range, geometric_mean)` when averaging over a time
  range).


## New extension mechanism allowing extension chains

- Many functions are now S3 generics with methods for
  `MizerParams` or `MizerSim` objects, and users can define their own subclass
  methods to modify mizer behaviour (#330).

- New composable extension chain infrastructure: `registerExtensions()`,
  `getRegisteredExtensions()`, `coerceToExtensionClass()`,
  `clearExtensionChain()`, and `registerExtension()`. Extension classes are S3
  marker classes; `MizerSim` derives its extension chain from
  `sim@params@extensions`. Extensions that do not provide a marker class remain
  metadata-only and do not trigger the S3 projection-rate dispatch path.

- S3 projection hooks have been added for all standard mizer rate functions.
  Extension-aware projections dispatch through `projectRates()`,
  `projectEncounter()`, `projectFeedingLevel()`, `projectEReproAndGrowth()`,
  `projectERepro()`, `projectEGrowth()`, `projectDiffusion()`,
  `projectPredRate()`, `projectPredMort()`, `projectFMort()`, `projectMort()`,
  `projectRDI()`, `projectRDD()`, and `projectResourceMort()` — while models
  without extensions continue to use the pre-resolved `mizerRates()` pipeline
  directly, with no per-step overhead.
  
- The `MizerSim` accessors `getParams()`, 
  `validSim()`, `N()`, `NResource()`, `finalN()`, `finalNResource()`,
  `idxFinalT()`, `getTimes()`, `getEffort()`, and are now
  registered as S3 generics with `MizerSim` methods, making extension-specific
  methods possible. `validParams()` is also now an S3 generic.

- `saveParams()` now serialises extension objects as plain `MizerParams`
  objects while preserving their extension chain, and `readParams()` restores
  the appropriate extension class. New `saveSim()` and `readSim()` helpers
  provide the same lifecycle for `MizerSim` objects.

- Extension installation support now integrates `pak` for managing missing or
  outdated extension packages.

- New vignette
  [Extending mizer](https://sizespectrum.org/mizer/articles/extending-mizer.html)
  documents when to use `setRateFunction()`, `setComponent()`, and
  `customFunction()`, summarises required function signatures and return shapes,
  and gives worked examples for both a custom encounter function and an added
  ecosystem component. A companion vignette
  [Using extension packages](https://sizespectrum.org/mizer/articles/using-extension-packages.html)
  is aimed at users of extension packages, and
  [Creating a mizer extension package](https://sizespectrum.org/mizer/articles/creating-extension-packages.html)
  guides extension authors through setting up a new extension package.

- `setRateFunction()` now validates the registered function by calling it with
  test inputs and checking that the return value has the correct dimensions,
  catching mismatched custom rate functions at registration time rather than
  during a simulation run. Closes #167.

- `setComponent()` now accepts optional `colour` and `linetype` arguments and
  applies them via `setColours()` and `setLinetypes()` so added components can
  be styled directly in plots.

- The `plot()` and `summary()` methods for `MizerParams`, `MizerSim`, and the
  mizer array classes are now registered as S3 methods rather than S4 methods,
  so `plot()` and `summary()` remain plain S3 generics when mizer is loaded,
  avoiding interference with S4 method dispatch for other packages.

## Species parameters for external mortality, encounter and diffusion rates

See the [model description](https://sizespectrum.org/mizer/articles/model_description.html) vignette for
the mathematical details.

- New species parameters `z_ext` (default 0) and `d` (default `n - 1`) add an
  optional power-law term to the external mortality: `mu_ext(w) = z0 + z_ext *
  w^d`. When `z_ext` is zero (the default) the behaviour is unchanged. Closes
  #329.

- New species parameter `E_ext` (default 0) sets the coefficient of the
  external encounter rate power law. `setExtEncounter()` now calculates the
  default external encounter rate as `E_ext * w^n` when no custom array is
  supplied, matching the pattern of `setMaxIntakeRate()`. A `reset` argument is
  also added to `setExtEncounter()` to force recalculation from species
  parameters.

- New species parameter `D_ext` (default 0) sets the coefficient of the
  external diffusion rate power law. `setExtDiffusion()` calculates the default
  array from species parameters when no custom array is supplied.

## Other improvements

- The `MizerSim` methods of the rate-getter functions (`getEncounter()`,
  `getFeedingLevel()`, `getEReproAndGrowth()`, `getERepro()`, `getEGrowth()`,
  `getDiffusion()`, `getPredRate()`, `getPredMort()`, `getMort()`, `getFMort()`,
  `getFMortGear()`, `getRDI()`, `getRDD()` and `getFlux()`) are now much faster.
  They resolve the rate functions and validate the parameters once and then, at
  each saved time step, calculate only the rates needed (and their
  dependencies) rather than re-resolving and recomputing the whole rate chain.
  The speed-up grows with the depth of the rate chain, e.g. roughly 100× for
  `getRDI()` and `getFlux()` on a 50-step simulation.

- New `scaleRates(params, factor)` function that rescales all rates in a model
  by a given factor. This is equivalent to a time rescaling: it speeds up or
  slows down all dynamics without affecting the steady state. All rate slots
  (`search_vol`, `intake_max`, `metab`, `mu_b`, `ext_encounter`,
  `ext_diffusion`, `catchability`, `rr_pp`) and their associated species
  parameters (`gamma`, `h`, `ks`, `k`, `z0`, `z_ext`, `z0pre`, `E_ext`,
  `D_ext`, `R_max`) are rescaled consistently.

- New `getTrophicLevel()` function returns a matrix (species × size) with the
  trophic level of individuals at each size, accounting for ontogenetic diet
  shifts by integrating the consumption-weighted average prey trophic level
  over the individual's growth trajectory. New `getTrophicLevelBySpecies()`
  returns the consumption-rate-weighted mean trophic level per species. Both
  functions accept `MizerParams` and `MizerSim` objects. Closes #307.

- New `expandSizeGrid()` function (an S3 generic) expands the size grid of a
  `MizerParams` object to a new minimum and/or maximum size while preserving
  all existing species data. Both `addSpecies()` and `expandSizeGrid()` now
  preserve the `MizerParams` subclass. `upgradeParams()` also preserves
  `MizerParams` subclasses and their extra slots.

- `compareParams()` output is now printed in a human-readable format, with each
  difference as its own block separated by blank lines. When array slots differ,
  the max absolute difference is shown per species. When slots differ only in
  their `comment` attributes, both comments are displayed. Closes #205.

- `summary()` for `MizerParams` and `MizerSim` now displays metadata from the
  `@metadata` slot, including title, description, authors, DOI, URL, mizer
  version, and creation/modification timestamps (when set). Closes #294.

- New `str()` methods for `MizerParams` and `MizerSim` objects, and the mizer
  array classes (`ArraySpeciesBySize`, `ArrayTimeBySpecies`, and
  `ArrayTimeBySpeciesBySize`), showing a clean, compact overview of their
  structures without dumping large amounts of internal data.

- A new `steady` argument to `addSpecies()` controls whether `steady()` is
  called after adding the new species.

- `constantEggRDI()` now accounts for diffusion across the egg-size boundary,
  including when `project()` uses the `"predictor-corrector"` method.

- `setRateFunction()` now validates custom RDI functions with the same
  `diffusion` argument that they receive during projection.

- Growth is now forced to always be non-negative, preventing unphysical
  shrinkage. No warning is issued when growth stops at or after maturity size.

- Added `info_level` argument to `projectToSteady()`, `steady()`, `setParams()`,
  `newCommunityParams()`, `newTraitParams()`, `matchBiomasses()`,
  `matchNumbers()`, `matchYields()` and `addSpecies()`to control the
  verbosity of information messages, consistent with `newMultispeciesParams()`.
  Set `info_level = 0` to suppress all messages. Closes #290. 

- `t_max` and `t_save` arguments in `project()` are now respected even when an
  effort array is supplied. When `t_max` is provided, the simulation extends
  beyond the times in the effort array using the last known effort values. When
  `t_save` is provided, it controls the save frequency with effort values
  interpolated as needed (#231).

- `getBiomass()` now has a `use_cutoff` argument to restrict the biomass
  calculation to sizes above the `biomass_cutoff` species parameter.
  `plotBiomass()` and `plotlyBiomass()` also gain this argument.

- `setResource()` now allows `resource_level = 1`. When balancing would
  otherwise divide by zero because the resource capacity equals the current
  resource abundance at positive consumption, the capacity is increased
  slightly with a warning instead of failing early.
  
- `project()` now warns when `t_max` is not a multiple of `t_save` and ensures
  that the state at `t_max` is always saved, even if the final save interval is
  shorter than `t_save`. (#341)
  
- New function `psi()` returns an `ArraySpeciesBySize` with the population-level
  reproductive proportion.

- `age_mat_vB()` is now exported.

- New [Cheatsheet: Analysis and Plotting](https://sizespectrum.org/mizer/articles/cheatsheet-analysis-and-plotting.html)
  vignette provides a quick reference for all functions that access simulation
  arrays, compute summaries, calculate indicators, and create plots.
  Closes #176.

## Bug fixes

- `getFMort()` on a `MizerSim` object was silently dropping the component
  names from `n_other` when passing it to the rate function and its
  dependencies (`getEGrowth()`, `getPredMort()`), causing failures whenever
  rate functions accessed `n_other` by name (e.g. `n_other[["resource"]]`).
  The implementation has been refactored to use the same `plyr::aaply` pattern
  as `getFeedingLevel()` and `getPredMort()`.

- `getFMort.MizerSim()` was not passing the time argument `t` to user-defined
  fishing mortality functions.

- `plotSpectra()` was incorrectly forcing the y-axis lower limit to 1e-20
  (instead of auto-scaling to the data) and was using `min(params@w) / 100`
  as the default lower w-axis limit even when `resource = FALSE`, where
  `min(params@w)` is more appropriate.

- `upgradeParams()` was silently dropping some slots (e.g. `resource_dynamics`)
  and was not preserving `MizerParams` subclasses and their extra slots when
  upgrading older objects.

- `getMeanMaxWeight()` now correctly applies the species selector to the
  denominator.

- `plotDataFrame()` now correctly applies custom log-scale x breaks.

- `get_size_range_array()` no longer gives an error when no size brackets are
  selected.

## Breaking changes

- The default `ratio` argument in `plotBiomassObservedVsModel()` and
  `plotlyBiomassObservedVsModel()` is now consistently `FALSE` for all object
  types. Calls that relied on the previous default ratio plot should now set
  `ratio = TRUE`.

- The first argument of `plotBiomass()`, `plotYield()`, `plotYieldGear()` and
  their `MizerSim` methods and `plotly*` wrappers has been renamed from `sim`
  to `object` for consistency with other plot generics. Calls using
  `sim = ...` as a named argument must be updated to `object = ...`.

- The names of the dimnames of the arrays returned by `getMort()`, 
  `getPredRate()` are now `sp` and `w` to be in line with other
  functions like `getFMort()`.
  
- Functions that return arrays of the form (species x size), (time x species)
  or (time x species x size) now return them with extra attributes and an S3
  class of `ArraySpeciesBySize`, `ArrayTimeBySpecies` or 
  `ArrayTimeBySpeciesBySize`. While this does not change their old behaviour,
  the differences will be flagged by functions like `is.identical()`.

- Because `plotDataFrame()` now correctly applies custom log-scale x breaks,
  the axis ticks in plots that use this function have changed.

- `plotDiet()` no longer accepts a `time_range` argument.

# mizer 2.5.4

- New function `renameGear()` to rename gears in a MizerParams object, similar 
  to `renameSpecies()`.
- `addSpecies()` now proceeds with a warning instead of an error when species
  growth stops after maturity (#315).
- `matchBiomasses()` and `matchNumbers()` now provide more informative error
  messages.
- `plotDiet()` now restricts the plot to size ranges with meaningful biomass
  density (#317).
- The `wlim` and `ylim` arguments in plotting functions now set the actual axis
  limits instead of just zooming (#320).
- The legend in `plotlyFeedingLevel()` is improved when critical feeding level
  is included.
- `species` and `gears` columns are now never factors, so no longer need to
  call `as.character()` so often.
- `validParams()` also calls `validGearParams()`.
- `validParams()` checks that `w_min` is valid for all species and increases it
  if necessary.
- `validSpeciesParams()` now also sets default for `p` to be equal to `n`.
- `species_params<-()` and `given_species_params<-()` now check that species
  names match.
- The `params` argument in `l2w()` and `w2l()` has been renamed to `species_params`
  to follow mizer's convention that `params` refers to a MizerParams object.

## Bug fixes

- `animateSpectra()` now uses consistent colours and preserves colour identity
  across frames (#321).
- `getReproductionProportion()` no longer returns incorrect proportions > 1 (#299)
- `setResource()` now correctly applies the `w_pp_cutoff` parameter to the 
  carrying capacity and initial resource abundance when changed without 
  providing `resource_capacity`(#306).
- Predation kernels are now truncated as documented.
- `given_species_params()` no longer makes unwanted changes to the species
  parameters.
- `steadySingleSpecies()` no longer changes `time_modified`.

# mizer 2.5.3

A patch update so that users who had changed `w_max` manually in their model
will not get unhelpful error messages when trying to use their model in the
new version. General checking of parameters is made more robust. In particular

- `validSpeciesParams()` has extra checks on consistency of species parameters
- `validParams()` checks that rate arrays contain finite numeric values
- `validSim()` checks that simulation results are finite and truncates the
  simulation if they are not.

# mizer 2.5.2

- Fixed bug that had led `newCommunityParams()` to set up resource parameters
  differently since version 2.4.0 (#293)
- `addSpecies()` now correctly preserves all `species_params` of the existing
  model.
- `addSpecies()` no longer requires new species to grow to maximum size, only
  maturity size is required.
- Now `validGivenSpeciesParams()` validates the given species parameters without
  adding defaults and `validSpeciesParams()` validates and returns a completed
  species parameter dataframe.
- New species parameter `w_repro_max` giving the size at which a species 
  invests 100% of its energy into reproduction. Set to `w_max` by default.
- `removeSpecies()` now also removes species parameters that are not set for
  any of the remaining species.
- Changing `w_max` now also correctly updates `ft_mask` (#296).
- `compareParams()` now also spells out differences in given species parameters.
- `getDiet()` now also includes the contribution of the external encounter rate
  to the diet.
- `setPredKernel()` now throws an error if some of the required predation kernel
  parameters are NA.
- In `plotYieldGear()` one can select a subset of gears with new `gears` 
  argument.
- New helper function `valid_gears_arg()` to check the `gears` argument in 
  functions that take a `gears` argument.
- Improved scaling of the y-axis in `plotGrowthCurves()`.
- `steadySingleSpecies()` no longer requires species to grow to `w_max`.
- `matchGrowth()` now also rescales the external encounter rate.
- `setExtEncounter()` no longer resets the external encounter rate to zero when
  called without the `ext_encounter` argument.
- The function `plotBiomassObservedVsModel()` now plots the ratio of modelled
  to observed biomass as default (`ratio = T`), as this is more useful visually
  to see how far off modelled biomass is from observed biomass.
- The `time_modified` field is now updated correctly by `steadySingleSpecies()`,
  `setColours()` and `setLinetypes()`.
- Deprecated `matchYields()` and `calibrateYield()`.
- Improved some unit tests.
- Some improvements to documentation.

# mizer 2.5.1

This is a patch release made necessary by a change in CRAN's requirement
regarding the vignettes. It also includes a bug fix:

- `project()` and `projectToSteady(..., return_sim = TRUE)` now correctly 
  returns also the other components of the MizerSim object stored in `n_other`.
  #285


# mizer 2.5.0

This release introduces a change that requires you to upgrade your old 
MizerParams and MizerSim objects with `upgradeParams()` or `upgradeSim()`.

## External encounter rate

Now the model can include an external encounter rate that represents the
rate at which a predator encounters food that is not explicitly modelled.
This encounter rate is set with `setExtEncounter()` or `ext_encounter<-()`
and can be read with `getExtEncounter()` or `ext_encounter()`. So this is
similar to how external mortality is handled.

## Given versus calculated species parameters

You can now use `given_species_params()` to see the species parameter
values that you have explicitly specified and `calculated_species_params()`
to see the species parameter values that mizer has calculated automatically or
set to defaults. You can continue to use `species_params()` to get all
species parameters, irrespective of whether they were given or calculated.

You can still set parameter values with `species_params<-()`, but you can also
use the stronger `given_species_params<-()` which not only sets the values you
give but also triggers a re-calculation of the calculated species parameters.
Using `given_species_params<-()` is therefore usually the better option.

## New mizer course

There is now a three-part mizer course at https://mizer.course.sizespectrum.org
with each part consisting of several tutorials, including code and exercises:

-   **Part 1: Understand**\
    You will gain an understanding of size spectra and their dynamics by exploring simple example systems hands-on with mizer.

-   **Part 2: Build**\
    You will build your own multi-species mizer model for the Celtic sea, following our example. You can also create a model for your own area of interest.

-   **Part 3: Use**\
    You will explore the effects of changes in fishing and changes in resource dynamics on the fish community and the fisheries yield. You will run your own model scenarios.


## Other improvements

- Warnings are given if user gives irrelevant species parameter values.
- Some messages have been converted to warnings and some to signals that are not
  shown as frequently.
- Frequent warnings are avoided when length-based and weight-based parameters 
  are both given and are inconsistent. #277
- Documentation of `effort` argument in `project()` is improved.
- An error message is given if a predation kernel returns negative values or
  is everywhere zero. #283

## Bug fixes

- When the coefficient `h` of the maximum intake rate is not given, it is now
  again given a default value. #282
- `matchGrowth()` no longer gives an error when there is no `w_inf` column. #279


# mizer 2.4.1

This minor release was made necessary to keep mizer on CRAN after a unit test
failed on macOS 13.3 with version 14.3 of the CLT toolchain.

# mizer 2.4.0

This release introduces a change that requires you to upgrade your old 
MizerParams and MizerSim objects with `upgradeParams()` or `upgradeSim()`.

See [mizer 2.4.0 blog post](https://blog.mizer.sizespectrum.org/posts/2022-12-23-mizer-240/)

## Avoid confusion between maximum size and von Bertalanffy asymptotic size

For an explanation see blog post at
https://blog.mizer.sizespectrum.org/posts/2022-11-30-dont-use-von-bertalanffy-growth-parameters/

The species parameter that specifies the size at which also the largest fish stop
growing is renamed from `w_inf` to `w_max`. The parameter `w_inf` is now 
reserved for the von Bertalanffy asymptotic size parameter. If you upgrade
your existing MizerParams object with `upgradeParams()` the `w_inf` column is
copied over to the `w_max` column automatically, but you may want to change
the values yourself if they do not currently reflect the maximum size of the
species. Otherwise the size distributions predicted by mizer will not match
observations.

## Set resource abundance rather than carrying capacity

The resource parameters `kappa` and `lambda` are now used to set the abundance
of the resource in the steady state rather than the carrying capacity, because
the latter is not observable.

While tuning the steady state using the `steady()` function the resource
abundance is now being kept fixed at the chosen value. Then the resource
dynamics can be switched on later with `setResource()` without changing the
steady state. At that stage you only choose either the resource intrinsic
growth rate or the resource carrying capacity and the other is determined by
`setResource()` in such a way that the resource replenishes at the same rate at 
which it is consumed. If you want to keep the old behaviour and switch off this
automatic balancing you have to add the `balance = FALSE` argument when calling
`setResource()`.

You can also choose between semichemostat dynamics `resource_semichemostat()`
or logistic dynamics `resource_logistic()` or you can write your own function 
implementing more sophisticated resource dynamics.

The `setParams()` function no longer includes the arguments for setting the
resource parameters. Instead you set these separately with `setResource()`.

## Automatically match growth rates

As explained in the blog post at https://blog.mizer.sizespectrum.org/posts/2022-11-30-dont-use-von-bertalanffy-growth-parameters/, 
the von Bertalanffy curves fitted to size-at-age
data are not suitable for estimating the size-dependent growth rates in mizer.
It is therefore now recommended that instead of von Bertalanffy parameters you
supply the age at maturity in the `age_mat` column of the species parameter
data frame. This is then used by mizer to calculate a default for the 
maximum intake rate parameter `h` if you do not supply this.

In the past, whenever you changed any model parameters, you needed to re-tune
other parameters to keep the growth rates in line with observations. There is
now a new function `matchGrowth()` that automatically scales the search volume,
the maximum consumption rate and the metabolic rate all by the same factor in
order to achieve a growth rate that allows individuals to reach their maturity
size by their maturity age while keeping the feeding level and the critical
feeding level unchanged. This function does not however preserve the steady
state, so you will need to also call `steady()` after matching the growth rates.


## Other improvements

* New function `steadySingleSpecies()` that only balances the size-spectrum
  dynamics while ignoring multi-species effects. In other words, it calculates
  the steady-state size spectrum of each species as it would be if the abundance
  of prey and predators could be kept constant at their current values.
* `plotGrowthCurves()` can now superimpose a scatterplot of size-at-age data
  if you supply this via the new `size_at_age` argument.
* New functions `calibrateNumber()` and `matchNumbers()` that are like
  `calibrateBiomass()` and `matchBiomasses()` but work with observed numbers
  instead of observed biomasses.
* New function `age_mat()` to calculate the age at maturity from the growth
  rate and the size at maturity.
* If an effort vector or effort array contains NA's, these are now replaced by
  the default effort value. #230
* The entries of the interaction matrix and of interaction_resource are no
  longer restricted to be less or equal to 1. #232
* If user supplies no row names in the interaction matrix but gives column names
  then the column names are also used as row names. #247
* `project()` now also works when called with a MizerSim object with additional
  components.
* `steady()` now preserves the RDD function in the MizerParams object rather
  than always setting it to "BevertonHoltRDD".
* When averaging abundances over time in `plotSpectra()` or `setInitialValues()`
  the user can now choose geometric averaging with `geometric_mean = TRUE`.
* The `w_mat25` species parameter is no longer filled in automatically if it is
  not supplied. This makes it easier to change `w_mat` without having to change
  `w_mat25` at the same time.
* `compareParams()` now also checks the validity of its second argument.
* Hide the printing of messages about chosen defaults in `newTraitParams()`.
* Higher values for the `info_level` argument in `newMultispeciesParams()` now
  leads to more messages.
* Giving more helpful messages in `validSpeciesParams()`. #136
* New helper functions `l2w()` and `w2l()` for converting between length-based
  and weight-based species parameters. #258
* Check that assessor functions for MizerSim slots are called with a MizerSim
  object.
* Add `style` argument to `plotDataFrame()` to facilitate producing area plots.
* Add `wrap_scale` argument to `plotDataFrame()` to control scaling of axes in
  faceted plots.
* `plotDiet()` can now show diets of several predator species in a faceted
  plot. #267
* Change from `size` to `linewidth` aesthetic to avoid warnings in new version
  of ggplot2.
* Better error message when functions are called with no valid species selected.
  #251
* If there are no differences then `compareParams()` says so clearly.
* `getReproductionLevel()` works as long as `R_max` is set. #252
* Converted several unit tests to edition 3 of testthat package.
* Improved documentation for `gear_params()`.
* Improved defaults can now be implemented while keeping backwards compatibility
  via `defaults_edition()`. #186
* New defaults edition 2: catchability = 0.3 instead of 1, initial effort = 1
  instead of 0. #243
* In defaults edition 2, `get_gamma_default()` ensures a feeding level of `f0`
  for larvae also if `interaction_resource` is not equal to 1. #238
* Set default linecolour and linetype for external mortality.

  
## Bug fixes

* Restored the line colours to `NS_params`
* Comment field now preserved by `set_species_default()`. #268
* Comment on `w_inf` no longer leads to error in `plyr::aaply()`. #269
* Can now again set `url` field in metadata. 
* Correct species now listed in the legend of `plotYieldObservedVsModel()` and
  `plotBiomassObservedVsModel()`. #266
* Standard order for legend in `plotDiet()` restored after change to `ggplot2`
  package. #265
* Fix handling of column names when interaction matrix is read from .csv file.
  #263


# mizer 2.3.1

* Resolved conflict in `mizerPredRate()` between the argument `t` and the 
  function `base::t()`.
* Assert that upgradeParams() must be called with a MizerParams object and 
  `upgradeSim()` with a MizerSim object.
* Errors changed to warnings in `getRequiredRDD()`
* `renameSpecies()` no longer fails when linecolour and linetype are of
  different lengths.
* matchYields() now also works for a model with only a single species.
* `setInitialValues()` can now average over a time_range.
* `getSSB()`, `getBiomass()`, `getN()`, `getYieldGear()` and `getYield()`
  can now be called with a MizerParams object as well as with a MizerSim 
  object. (#200)
* Updated the shiny app in inst/shiny/selectivity_effects to current mizer
  version.

# mizer 2.3.0

## New features

* New plots `plotBiomassObservedVsModel()` and `plotYieldObservedVsModel()`
  contributed by @SamikDatta., together with their plotly counterparts.
* New `calibrateBiomass()`, `calibrateYield()` to set the model scale to agree
  with total observed biomass or total observed yield. Uses the new
  `scaleModel()`.
* New `matchBiomasses()` and `matchYields()` will try to adjust the abundances
  of the species to produce the observed biomasses or yields.
  See blog post at https://blog.mizer.sizespectrum.org/posts/2021-08-20-a-5-step-recipe-for-tuning-the-model-steady-state/ .
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
* mizer can be [run on binder](https://mybinder.org/v2/gh/sizespectrum/mizer/HEAD?urlpath=rstudio)

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
* Mizer re-exports the `melt()` function from the reshape2 package which allows
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
