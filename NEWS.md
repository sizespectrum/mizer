# mizer 3.2.1.9000

This development version adds experimental tools for analysing the dynamic
stability of steady states, tools for checking whether a model actually is at
its steady state, a general `sizeIntegral()` for integrals over the size
spectrum, and a substantial extension of the plotting functions, in particular
of plots against a length axis. It also routes nearly everything mizer says
while building or changing a model through a single mechanism controlled by
`info_level`.

## Dynamic stability

- New experimental `steadyNewton()` finds a steady state by solving the
  steady-state equation directly with a Newton-type root finder (via the
  `nleqslv` package) instead of running the dynamics to convergence. Unlike
  `steady()` it converges even when the steady state is dynamically unstable, and
  it discovers the support of the steady state automatically. Currently supports
  the default semichemostat resource dynamics.

- New experimental `getStability()` analyses the dynamic stability of a mizer
  steady state by computing the eigenvalues of the linearised one-step-ahead map
  at the fixed point. To eliminate the artificial temporal numerical diffusion
  introduced by the implicit solver, it maps the discrete numerical eigenvalues
  back to their exact continuous-time equivalents via $\lambda = (1 - 1/\mu) / \Delta t$.
  It reports whether the steady state is stable or unstable (based on the real
  parts of the continuous eigenvalues), the maximum real part, the `spectral_radius`
  of the discrete numerical solver for a given time step `dt`, and — when the system
  emergent limit cycle. By default the resource is treated as a
  quasi-static fast variable (valid for semichemostat dynamics); setting
  `include_resource = TRUE` gives the full coupled (fish + resource) Jacobian,
  useful for verifying that the quasi-static approximation makes little difference.
  The stability list also includes
  `leading_eigenvectors`: a complex array `(n_species, n_sizes, 2)` of the
  top two eigenvectors reshaped into the fish abundance space, normalised to
  maximum modulus 1. The rate functions are only ever evaluated at states
  satisfying `N >= 0` — where a centred difference would push a cell negative,
  the column is differenced forwards from the unperturbed state instead — so a
  custom rate function registered with `setRateFunction()` never has to be
  defined at negative abundances.

- New experimental `getLimitCycleSim()` takes the output of `steadyNewton()`
  and constructs a `MizerSim` covering one period of the limit cycle in the
  linear approximation. The trajectory is
  \eqn{N(t) = N^* + A\,\text{Re}[e^{i\theta t}\,\mathbf{v}]}, where
  \eqn{\mathbf{v}} is the leading complex eigenvector and the amplitude \eqn{A}
  is scaled so the maximum relative perturbation equals the `amplitude`
  argument (default 10\%). The returned object can be passed directly to
  `plotBiomass()`, `plotSpectra()`, and other standard mizer plot functions.

- New experimental `plotBifurcation()` draws a bifurcation diagram over fishing
  effort. For each effort value it follows the attractor of the full dynamics
  and plots the long-term range of a summary quantity (biomass, yield or SSB).
  A stable steady state appears as a single line and a limit cycle as a band
  between the minimum and maximum, so a Hopf bifurcation shows up as the effort
  at which the band opens up. The settling stage runs `projectToSteady()`, whose
  `tol`, `amplitude_tol` and `extinction_threshold` are exposed for tuning.

- `steady()` and `projectToSteady()` now report the nature of the solution they
  converged to via a `"convergence"` attribute on the returned object (mirroring
  the `"stability"` attribute of `steadyNewton()`). It records whether the run
  settled on a stable steady state, a limit cycle, or neither, together with the
  cycle period and relative amplitude when a cycle is found. Limit cycles are
  detected from a per-species biomass series sampled at the new `t_save`
  resolution (default `dt`), so detection no longer relies on the cycle period
  being commensurate with `t_per`. The relative-amplitude floor for calling an
  oscillation a limit cycle is a separate `amplitude_tol` argument (default
  `0.01`), independent of the fixed-point convergence tolerance `tol`, and a
  species is treated as extinct once its reproduction falls below the
  `extinction_threshold` fraction (default `1e-6`) of its value at the start of
  the run.

- `getStability()` and `getLimitCycleSim()` warn when they are handed a model
  that is not at a steady state. Both linearise the dynamics *at* the stored
  state, so on a state that is not a fixed point their eigenvalues describe the
  neighbourhood of a point the model is not sitting at (#495).

## Steady state and calibration

- New experimental S3 generic `isSteady()` returns `TRUE` if a model is at its
  steady state (within a tolerance, defaulting to 0.05/year), `FALSE` otherwise.
  Extension packages can provide custom methods.

- New experimental `getSteadyResidual()` answers the question every calibration
  workflow otherwise has to remember to ask: *is this model still at its steady
  state?* It returns the rate at which each species' abundance would change if
  the model were projected forward, as a per-capita rate in 1/year, so zero
  means the model is on a fixed point. For the consumers the value is exact
  rather than a finite difference: the backward-Euler transport coefficients
  used by `project()` satisfy `A N - S = -dt dN/dt` identically. Everything is
  evaluated with the model's own reproduction function and its own
  `resource_dynamics`, so unlike `steadyNewton()` it works whatever the resource
  dynamics are. The result is an `ArraySpeciesBySize`, so
  `plot(getSteadyResidual(params))` shows which species and which sizes have
  moved (#495).

- `summary()` of a `MizerParams` object now reports the model's biomass drift
  and whether that counts as being at a steady state. Whether a model is settled
  is not visible from any of the other parameters shown, and getting it wrong is
  the most common way a calibration goes quietly wrong (#495).

- `project()` gains an experimental `check_steady` argument. With
  `check_steady = TRUE` it warns when it is handed a model that is not at its
  steady state, which catches the mistake of forgetting to re-run `steady()`
  after a `match…()` step. It defaults to `FALSE`, because projecting a model
  away from its steady state is a perfectly normal thing to do. The check is
  made at the effort stored in the params object rather than at the effort
  passed to `project()`, so running a fishing scenario at a new effort does not
  warn (#495).

- The `"convergence"` attribute attached by `projectToSteady()` and `steady()`
  gains a `residual` field giving how far the state reached actually is from a
  fixed point. The `distance` field only compares two states `t_per` apart on
  whatever scale the distance function uses; `residual` measures the thing
  itself, and `steady()` now says when the two disagree — a run declared
  converged whose biomasses are still visibly moving (#495).

- `matchBiomasses()`, `matchNumbers()` and `matchGrowth()` now say that they
  have moved the model off its steady state, turning an
  instruction the documentation had to keep repeating into something the package
  says at the moment it becomes true. The `calibrate…()` functions and
  `scaleModel()` deliberately do not: they apply one overall scaling factor,
  which is an exact symmetry of the model and leaves the steady state untouched
  (#495).

- New `reproduction_level()` and `reproduction_level<-` accessor and replacement
  functions allow reading and changing the reproduction level while preserving
  the steady state, matching the syntax of `resource_level()` and
  `resource_level<-`. `getReproductionLevel()` is deprecated in favour of
  `reproduction_level()`.

- `matchGrowth()` gains the `info_level` argument that the other `match…()`
  functions already had, and `matchBiomasses()` and `matchNumbers()` now
  actually honour theirs.

- The biomass and the number variants of the calibration and matching functions
  now share one implementation each. `calibrateBiomass()` and
  `calibrateNumber()` differ only in whether the size integral carries a factor
  of the weight, as do `matchBiomasses()` and `matchNumbers()`, but each was a
  separate copy of the code, and the copies had drifted apart. The four
  functions behave exactly as before on the default quadrature scheme; what
  changes is that a correction to how an observation is compared with the model
  is now made once instead of once per variant (#504).

## Integrals over the size spectrum

- New experimental `sizeIntegral()` calculates any integral
  \eqn{\int N_i(w)K_i(w)dw} over the size spectrum. It is now the recommended
  way to write your own summary or indicator function: it selects the size
  range, applies the quadrature scheme the model is actually on and wraps the
  result in the appropriate mizer array class, so none of those rules need to
  be remembered. It takes the weighting factor \eqn{K} via the `weighting`
  argument in any of the shapes mizer's own arrays come in, from a single
  number to a gear x species x size array, and keeps the extra dimensions in the
  result. `getBiomass()`, `getN()`, `getSSB()`, `getYield()`, `getYieldGear()`
  and `getProportionOfLargeFish()` are now all implemented with it (#494).

- New experimental `bin_average_weight()` prepares the weight of an integral
  over the size spectrum so that the integral uses the quadrature scheme the
  model is actually on. Use it when writing your own indicator or diagnostic
  function: it is gated on the `bin_average` entry of `second_order_w()`, so an
  indicator written with it is unchanged on the default scheme and correct to
  second order when second-order bin-averaging is switched on. Previously only
  available internally as `bin_average_summary_weight()`.

- New experimental `encounter_kernel()` returns the predation kernel that
  `mizerEncounter()` actually uses under whichever quadrature scheme the model
  is on, to be paired with the plain point prey weight `w_full * dw_full`. Any
  diagnostic that decomposes the encounter rate must use this rather than
  `getPredKernel()`, which is point-sampled on the grid and is intended for
  plotting and for supplying a custom kernel.

## Plotting

- Mizer arrays now state what kind of quantity they hold. Every array
  constructor gains a `type` argument: `"value"` (the default) for a rate or an
  amount, `"density"` for an amount per gram of body weight, `"proportion"` for
  a fraction. Two things follow from it.

  A `"density"` is multiplied by the appropriate Jacobian when it is plotted
  against a length axis (`size_axis = "l"`), and its units are restated from
  `1/g` to `1/cm`. This replaces the guess mizer used to make from the array's
  name and units, which recognised only densities that happened to be called
  "Number density" or to have units "1/g" — and so missed `getFluxGradient()`,
  whose values were plotted per gram against a length axis and labelled
  `g^-1/year`. Arrays that declare no type still fall back to that guess, so
  existing code and saved objects are unaffected.

  A `"proportion"` — `getFeedingLevel()`, `getCriticalFeedingLevel()`,
  `maturity()`, `repro_prop()`, `psi()`, `resource_level()` — is plotted on a
  linear y axis showing the whole of the interval from 0 to 1, so the value can
  be read against the scale it belongs to. The range is only ever widened to
  include the data, never narrowed to that interval: the critical feeding level
  and the resource level can both legitimately exceed 1, and their plots show
  it. `plot(getFeedingLevel(params))` therefore now shows the same y range that
  `plotFeedingLevel()` always has.

- `plot()`, `plot2()`, `addPlot()` and `animate()` on an array that holds a
  density gain a `per_log_size` argument, which expresses the values per
  logarithmic size rather than per size. This is the same change of measure that
  `size_axis` makes — both rescale the density by a Jacobian — so the two now
  sit side by side, and `plotSpectra()` is no longer the only way to see a
  spectrum per log size. Unlike `size_axis` it needs no weight-length
  relationship, so the resource classes take it too. Asking for it on an array
  that does not hold a density is now an error; it used to be swallowed silently
  by `...`.

- `plotSpectra()`, `plotSpectra2()`, `plotCDF()`, `plotCDF2()` and `animate()`
  now let you choose the plotted quantity with two independent arguments
  instead of the single `power`: `biomass` selects a biomass density rather
  than a number density and the new `per_log_size` selects a density with
  respect to logarithmic size rather than with respect to size. The power of
  the weight is the sum of the two, which is why `power = 1` was ambiguous: it
  is both the biomass density and the number density in log size, and mizer had
  to assume the former when choosing the y-axis label and the Jacobian for a
  length axis. The number density in log size, `plotSpectra(sim, biomass =
  FALSE, per_log_size = TRUE)`, is therefore now available for the first time
  with the correct label and length conversion. The `power` argument keeps
  working as before and remains the only way to ask for a power that is not the
  sum of the two flags, but supplying it together with a contradictory flag is
  now an error instead of silently ignoring the flag (#501). `plotCDF()` does
  not accept `per_log_size`, because the cumulative distribution does not
  depend on that choice. As a side effect, `plotlySpectra()`, `plotlyCDF()`,
  `plotlySpectra2()` and `plotlyCDF2()` now honour `biomass`, which they used
  to drop because they always passed `power` on internally.

- The resource can now be shown on length-based plots. `resource_params()`
  gains the weight-length parameters `a` and `b`, defaulting to the equivalent
  spherical diameter of an organism with the density of water, `a = pi/6` and
  `b = 3`, which is the convention plankton ecology uses for a composite of many
  taxa. `plotSpectra(params, size_axis = "l")` therefore includes the resource
  spectrum, where it used to drop it silently, and the resource array plots and
  `animate()` gain `size_axis` and `llim`. The parameters feed none of the
  rates. Note that the resource then sits on the length axis at its own
  convention: a fish of a given weight is about 3.7 times longer than a sphere
  of that weight. That difference is real rather than an artefact — a 1 mg
  copepod really is shorter than a 1 mg fish larva — but it does mean the
  resource and the species are measured differently.

- The total line is now shown on length-based plots, where it used to be
  dropped. A total can only be formed once every line sits on the same
  coordinate, and on a length axis the lines do not: each species, and the
  resource, converts weight to length with its own allometric relationship, so
  at a given length they sit at different weights. The total is therefore
  summed *after* the conversion — at equal length rather than at equal weight —
  interpolating each series onto the union of all the size coordinates,
  logarithmically in size, with a series contributing nothing outside its own
  range. Where the series already share a grid, which is always the case on a
  weight axis, the union is that grid and nothing is approximated: the
  weight-axis total is unchanged.

  `plotSpectra2()`, `plotSpectraRelative()`, `plotCDF()` and `plotCDF2()` keep
  their total on a length axis too. The comparison plots used to convert the
  axis after assembling the two spectra, so the total they had been given —
  summed at equal weight — arrived at the conversion with no species to convert
  it by and was dropped. They now let `plotSpectra()` convert, so the total they
  receive is already the total on the axis being plotted. As a result
  `plotSpectra2()` also applies `ylim` the way `plotSpectra()` does on a length
  axis: it could not before, because the values it was filtering were a Jacobian
  away from the ones the limits described, so with `return_data = TRUE` it
  returned values outside the limits that a single spectrum plot would have
  dropped.

  `total = TRUE` also now means the same thing everywhere: the total of
  everything the object holds. For `plotSpectra()` that was already so — the
  resource and every species, whatever is drawn — and it stays so. The array
  plots have been brought into line: `plot(<array>, total = TRUE)` used to sum
  only the species that were selected for display, and now sums the whole
  array, so a plot of two species can be read against the community total.

- The array-plotting toolkit now covers every mizer array class. The resource
  classes `ArrayResourceBySize` (as returned by `getResourceMort()`,
  `resource_rate()`, `resource_capacity()`, `resource_level()` and
  `initialNResource()`) and `ArrayTimeByResourceBySize` (as returned by
  `NResource()` on a `MizerSim`) gain `plot2()`, `plotRelative()` and
  `addPlot()` methods, and `ArrayTimeByResourceBySize` also gains an `animate()`
  method, so you can now compare resource spectra before and after a model
  change or play one through a simulation the same way you already could for
  species rates. `ArrayTimeBySpeciesBySize` gains the `addPlot()` method it was
  missing. The `species`, `total` and `background` arguments do nothing for a
  resource array, which holds a single spectrum, so the resource methods warn if
  they are set (#468).

- `plotFeedingLevel()`, `plotlyFeedingLevel()` and `plotYield()` now delegate
  directly to the array plotting methods — `plot()` / `plotHover()` on
  `ArraySpeciesBySize` and `plot(getYield(object))` — so they share the
  behaviour and the arguments of those methods. `plotFeedingLevel()` keeps full
  support for `include_critical = TRUE` and its non-clipping proportion
  coordinate scaling.

- `plotYieldObservedVsModel()` gains a `gear` argument that restricts the
  comparison to the catch of the selected gears. Both the model yield and the
  observed yield are then taken from those gears only, so in a model where
  several gears catch the same species you can check the gears against their
  own observations instead of only their total. Without the argument the plot
  keeps comparing the yield summed over all gears. Because the species
  parameter `yield_observed` is a total over all the gears, per-gear
  observations have to be given in `gear_params()` (#286).

## Reporting information to the user

- New `mizer_info_level` option sets how much mizer tells you about the choices
  it makes, without your having to pass `info_level` to each call.
  `options(mizer_info_level = 0)` quietens mizer as a whole, including the
  functions that have no `info_level` argument of their own, such as
  `species_params<-()` and the rate setters. The `info_level` argument still
  overrides it for a single call, and its default is now
  `default_info_level()`, which reads the option.

- The "Creating a mizer extension package" article is now generated from a new
  `create-extension-package` skill, and is named `guide-create-extension-package`
  like the other guides; the old address redirects. Everything in the
  extending-mizer guide that only matters once you share an extension — method
  dispatch, bundled data objects, reporting to the user, upgrading saved objects
  — moved into it, leaving that guide to the mechanisms themselves. Its advice on
  marker classes has also been corrected: it still told you to define them with
  `setClass()`, which mizer 3.2 made unnecessary and which prevents your package
  from being chained with another.

- The reporting mechanism is now exported, so that an extension package can tell
  the user what it decided on their behalf through the same channel mizer uses,
  and have it obey the same switch:
  `signal_info()` raises a report, `with_info_level()` collects the reports
  raised inside a call and gives them together when it finishes,
  `signal_not_recalculated()` is the standard report for a setter that left a
  hand-set array alone, and `default_info_level()` reads the `mizer_info_level`
  option. Give your own constructors and setters an
  `info_level = default_info_level()` argument and forward it, rather than
  hard-coding a value in the call to `newMultispeciesParams()`, which would
  collide with a user's own `info_level`.

- The information mizer gives while it sets up or changes a model is now raised
  through one function, `signal_info()`, which says which quantity the report is
  about, how important it is, whether it is a message or a warning, and whether
  it should still be shown when nothing is collecting reports. The collecting
  handlers now nest by themselves, so a function can report the information
  raised inside it without having to know whether its caller is already doing
  so, and two different things said about the same quantity are both reported
  where previously the second overwrote the first. This is the machinery behind
  the frozen-rate warnings described below, and it is available to packages that
  extend mizer.

- Nearly every message and warning that mizer gives while building or changing
  a model now goes through that mechanism, including the ones in `steady()`,
  `steadyNewton()`, `projectToSteady()`, `validParams()`, `setInteraction()`,
  `setReproduction()`, `setBevertonHolt()`, `setResource()`, `newTraitParams()`,
  `newSingleSpeciesParams()`, `plotYieldObservedVsModel()` and the upgrade of
  an old object. They are collected into a single report rather than a stream,
  and `info_level` (or the `mizer_info_level` option) controls all of them
  alike, where before each function decided for itself. `steady()`,
  `projectToSteady()` and `validParams()` no longer implement
  their own `info_level` threshold, and `newSingleSpeciesParams()`,
  `setBevertonHolt()` and `steadyNewton()` gain an `info_level` argument.

## Species, gear and resource parameters

- Changing a species parameter that feeds a rate array you have set by hand now
  warns you that the change has no effect on the model. Previously
  `given_species_params<-()` recorded the new value in the species parameter
  table but left the model unchanged, and said nothing: the rate setters do
  emit a message in that situation, but the setter runs `suppressMessages()`
  over the recalculation to quieten the routine chatter, so the message never
  reached the user. The report is now a warning, which survives that, and it
  names the species parameters that were ignored, the quantity that is holding
  them back, and the call that puts the quantity back under the control of the
  species parameters, for example `setMetabolicRate(params, reset = TRUE)`. It
  is raised only when a parameter that actually feeds the frozen quantity
  changed, so models that mizer itself freezes arrays in, like those from
  `newTraitParams()` and `newCommunityParams()`, do not warn about unrelated
  changes (#489).

- Changing a resource parameter that feeds a resource array you have set by
  hand now warns you that the change has no effect, the same way a species
  parameter change does. `resource_params(params)$kappa <- ...` with a
  `resource_capacity()` that was set manually previously changed the stored
  scalar and nothing else, in complete silence. `setResource()` also now says
  when it leaves a frozen resource array alone although the resource
  parameters ask for a different value (#489).

- The division of labour between the two species parameter setters is now
  clear-cut: `given_species_params<-()` is the one that warns when a change you
  asked for cannot take effect, and `species_params<-()` stays quiet. That
  covers all three such diagnostics — a parameter overridden by another one you
  have given, a parameter feeding a rate array you have set by hand, and a gear
  parameter that mizer reads from `gear_params()`. Use `species_params<-()` in
  scripts and `given_species_params<-()` interactively, where the diagnostics
  are worth having (#496).

- Those three diagnostics now agree about what counts as a change. Clearing a
  given species parameter to `NA` is one: it hands the parameter back to
  mizer's calculation, and if the rate array that calculation feeds has been
  frozen, that instruction cannot be carried out and you are now told so.
  Adding a column that holds only `NA` is not a change, and no longer draws a
  warning that nothing else agreed with. Previously each of the three
  diagnostics answered the question differently, so an all-`NA` new column
  warned about a frozen array while clearing a value that was actually given
  was reported by none of them (#524).

- `given_species_params<-()` now warns when `k_vb` is specified on a model where
  `h` or `age_mat` is already given and will override `k_vb`.

- `setExtMort()` now warns when explicitly supplied `z0pre` or `z0exp`
  arguments are ignored because `z0` is already present in
  `given_species_params()` for every species or because `ext_mort` was
  supplied. A `z0` value that is present only in `species_params()` is now
  recognised as calculated and is recalculated from `z0pre`, `w_inf` and
  `z0exp`. When either argument is explicitly supplied and used, the resulting
  values are recorded in `given_species_params()`. Values calculated from the
  arguments' defaults remain calculated parameters and are not recorded there
  (#493).

- `setParams(params, reset = TRUE)` is now a documented argument of
  `setParams()` rather than something that happened to be forwarded through
  `...`. It thaws every rate array that `setParams()` sets.

- `get_gamma_default()` and `get_ks_default()` now provide informative error
  messages when called on a model with `h = Inf`, indicating that `gamma` or
  `ks` must be supplied explicitly.

- `gear_params<-()` and `calc_selectivity()` now list the missing column names
  in the error message when arguments needed for a selectivity function are
  missing from the gear parameter data frame.

## Other new functions

- New `knife_edge_length()` selectivity function that applies a knife-edge cut
  at a given **length** rather than a weight. Set `sel_func = "knife_edge_length"`
  and provide a `knife_edge_length` column in `gear_params()`. The length is
  converted to a cut-off weight via the length–weight parameters `a` and `b`.

- New `gaussian_mixture_pred_kernel()` supports multimodal feeding preferences
  represented by mixtures of Gaussian distributions on the log
  predator/prey mass-ratio scale.

## Other improvements

- `validParams()` is now about 15 times faster on an object that is already
  valid (#461). It used to redo its full work on every call, which made it
  expensive to apply the "validate at the boundary" principle consistently. The
  repair work — rebuilding the species parameter tables and the `w_min_idx` and
  `ft_mask` slots, and checking the structural validity of the object — is now
  skipped for an object that has already been through it. Mizer recognises such
  an object by a fingerprint calculated from the contents of the slots that the
  repair and the structural checks depend on. The fingerprint is recalculated on
  every call and is not stored on the object, so it cannot go stale: any change
  to any of those slots, made by any route including a direct slot assignment,
  gives a new fingerprint and triggers the full validation. The checks for
  non-finite values in the rate arrays are still made on every call, because
  they catch what the fingerprint cannot see. `validSim()` benefits
  automatically, because most of its cost was its nested `validParams()` call.
  One consequence: any warning or message that the repair issues about a
  condition it does not itself fix (for example that a species has a maximum
  size larger than the largest size in the model) is now issued only the first
  time an object with that content is validated in a session.

## Deprecations

- `matchYields()` and `calibrateYield()` have been removed. They were deprecated
  in mizer 2.6.0 and no use case for them was reported. Both worked by
  multiplying the abundance of a species at all sizes by a constant factor,
  which is the wrong lever for a yield: use
  `mizerExperimental::matchYield()`, which adjusts the catchability instead
  (#526).

- Eleven accessors that returned a rate array stored in the MizerParams object
  had two names that did exactly the same thing. The `get`-prefixed name is now
  soft-deprecated in favour of the bare name, which is the one that also has a
  replacement function (`catchability(params) <- value` and friends):
  `getCatchability()` → `catchability()`, `getSelectivity()` → `selectivity()`,
  `getInitialEffort()` → `initial_effort()`, `getPredKernel()` →
  `pred_kernel()`, `getSearchVolume()` → `search_vol()`, `getMaxIntakeRate()` →
  `intake_max()`, `getMetabolicRate()` → `metab()`, `getExtMort()` →
  `ext_mort()`, `getExtEncounter()` → `ext_encounter()`,
  `getMaturityProportion()` → `maturity()` and `getReproductionProportion()` →
  `repro_prop()`. The old names keep working; they warn once per session in
  code you run directly. The `get` prefix is now reserved for the functions
  that *calculate* something from the current state of a model, like
  `getEncounter()` or `getFMort()`.

- `getReproductionLevel()` is deprecated in favour of the new
  `reproduction_level()`.

- The `sim2` argument of `plotYield()` is deprecated in favour of
  `plot2(getYield(sim1), getYield(sim2))`.

## Bug fixes

- `steadyNewton()` no longer fails with a non-finite residual error when the
  initial guess contains exactly zero abundance in the active size range (as
  can happen when averaging a simulation over a limit cycle). It now scales the
  residual by a strictly positive initial guess instead.

- `projectToSteady()` limit-cycle detection now uses only the second half of the
  series (or the most recent 20 samples if the series is short) for its
  autocorrelation step, instead of the entire run. This stops large initial
  transients from drowning out a cycle that settles later.

- Setting `f0` to a value outside the interval `[0, 1)` now gives an immediate
  error, whether or not `gamma` has been supplied. Previously `f0 = 1`
  silently produced an infinite `gamma` and a non-finite `search_vol` when
  `gamma` was calculated, while an invalid `f0` supplied alongside `gamma`
  could be accepted and ignored (#517).

- The default values for the `gamma` and `f0` species parameters are no longer
  corrupted by a search volume that you have set by hand. `get_gamma_default()`
  measures the energy available to a predator by giving it a search volume
  coefficient of 1, but it obtained that search volume by calling
  `setSearchVolume()`, which refuses to recalculate a `search_vol` array you
  have set yourself. So mizer's own internal calculation was blocked along with
  yours and the available energy was read off your array instead, making the
  resulting `gamma` wrong by however much your array differed from the
  unit-gamma one — many orders of magnitude in realistic models. Both
  `get_gamma_default()` and `get_f0_default()` now compute the search volume
  they need directly from the species parameters, so they are unaffected by a
  frozen `search_vol` and remain exact inverses of each other (#488).

- `get_f0_default()` now respects `interaction_resource` under
  `defaults_edition() >= 2`, making it the exact inverse of
  `get_gamma_default()` also in a model where the species do not all feed on
  the resource with the same strength.

- Setting `f0` to a value outside the interval `[0, 1)` now gives an immediate
  error, whether or not `gamma` has been supplied. Previously `f0 = 1`
  silently produced an infinite `gamma` and a non-finite `search_vol` when
  `gamma` was calculated, while an invalid `f0` supplied alongside `gamma`
  could be accepted and ignored (#517).

- Changing the resource power-law parameters now refreshes the species search
  volume parameters that mizer calculated from them. Changing `lambda`
  recalculates calculated `q` and `gamma` entries, while changing `kappa`
  recalculates calculated `gamma`; explicitly supplied species values remain
  protected. Previously `resource_params(params)$lambda <- ...` and the
  corresponding `setResource()` calls rebuilt the resource arrays but silently
  left these calculated species parameters and `search_vol` stale (#497).

- `species_params<-()` no longer records a default that mizer filled in as a
  given species parameter. A parameter that the model does not yet carry and
  that `validSpeciesParams()` supplies a default for — most commonly the
  length-weight parameters `a` and `b`, which a model built without length data
  does not have — looked like a column the user had just added, and so was
  recorded as given and frozen against every later recalculation. On
  `NS_params` any assignment to `species_params()`, even one setting an
  unrelated parameter, added 24 given entries this way. Values you do supply,
  including in a column that is new to the model, are recorded as before
  (#496).

- `$` on a `species_params` or `gear_params` table no longer partially matches
  column names. In a model without length-weight parameters,
  `species_params(params)$a` returned the `alpha` column and `$b` the `beta`
  column, with the species names attached, so code converting weights to
  lengths silently got the assimilation efficiency and the preferred
  predator/prey mass ratio instead. Writing was never partially matched, so
  reads and writes disagreed about which column `$b` meant. A name that is not
  a column now gives `NULL`, with a warning naming the column that used to be
  returned (#487).

- `w_min` is now included in `given_species_params`, so the `min_w` argument to
  `newMultispeciesParams()` and `emptyParams()` is preserved across any operation
  that rebuilds the species parameters from the given ones. Previously, a
  `given_species_params<-` round-trip would silently reset `w_min` to 0.001 when
  `min_w` was smaller than that, and emit a spurious warning when `min_w` was
  larger (#460).

- `setParams()` now gives an error when it is passed an argument that none of
  the setter functions it calls accepts. Every one of those setters declares its
  `...` as unused, so any misspelled or misplaced argument was silently
  accepted and ignored: `setParams(params, metabolic = 99)` did nothing and said
  nothing. Arguments that belong to another function are named in the error
  together with the function that does accept them. In particular
  `setParams()` never called `setResource()`, so `setParams(params,
  resource_rate = ...)` had no effect on the model; the error now points at
  `setResource()`, and the deprecation warnings for `setResource(r_pp)` and
  `setResource(kappa)`, which used to recommend `setParams()`, now recommend
  `setResource()` too. `setResource()` likewise now errors on an argument it
  does not have, instead of ignoring it (#491).

- `info_level = 0` now really does silence all the information about default
  values, including the reports that until now were plain messages: the notes
  about the interaction matrix dimnames, about `no_w` being increased, about
  negative resource abundances, about what an upgrade changed, and the
  consistency corrections to `w_mat`, `w_mat25`, `w_min` and `w_repro_max`.
  Previously an information signal whose level was above `info_level` was
  passed on rather than dropped, so it could still be reported by a handler
  further out, and a plain message ignored `info_level` altogether.

- `repair_params()` now suggests `adjustSizeGrid()` instead of the deprecated
  `expandSizeGrid()` when `w_min` is smaller than the grid minimum.

- `compareParams()` now uses relative tolerance when comparing species
  parameters, so small-magnitude parameters like `gamma` (~1e-8) are no longer
  silently treated as equal when they differ by up to ~10%.

- `calibrateBiomass()`, `calibrateNumber()`, `matchNumbers()`,
  `plotBiomassObservedVsModel()` and `plotYieldObservedVsModel()` now integrate
  over the size grid with `sizeIntegral()` like everything else in mizer. They
  had each hand-rolled the sum, so they stayed on the first-order quadrature and
  cut the size range at a bin boundary even in a model with
  `second_order_w(params)["bin_average"]` switched on. In such a model a species
  matched to its observed biomass was then plotted off the 1:1 line, because the
  match and the plot measured the biomass differently, and the model yields in
  `plotYieldObservedVsModel()` came out 10-20% below `getYield()`, so the plot
  and its total-relative-error caption reported an under-prediction of the
  yields that was not there. Results on the default quadrature scheme are
  unchanged (#504, #529).

- `matchNumbers()` no longer re-tunes the reproduction parameters when it has
  nothing to match. Its guard against an empty selection of species never
  fired, so with no observations, or none for the species asked for, it left
  the abundances alone but still adjusted the reproduction parameters and
  reported a rescaling that had not happened (#504).

- `plotYieldObservedVsModel()` now finds the observed yield where mizer says it
  belongs. `yield_observed` is documented as a `gear_params()` column, and
  `given_species_params<-()` tells you to put it there, but the plot read it
  only from the species parameters, so a model that followed the advice got the
  error that no `yield_observed` had been provided. The plot now takes the
  observations from the gear parameters, summed over the gears catching each
  species, and still accepts them in the species parameters for the species
  that have no gear observation (#526).

- `getProportionOfLargeFish()` called on a `MizerParams` object gave a wrong
  answer. It multiplied the species x size abundance array by the vector of
  weights, which R recycles down the columns of the array instead of along the
  size axis, so every species but the first was weighted by the wrong sizes.
  The `MizerSim` method was unaffected, and the two now agree when applied to
  the same state (#494).

- `getN()` now applies the model's quadrature scheme to the size range it is
  given, so that under `second_order_w(params) <- c(bin_average = TRUE)` the
  bin straddling `min_w` or `max_w` contributes only partially, as it already
  did in `getBiomass()`. Results on the default first-order scheme are
  unchanged (#494).

- `getDiet(proportion = FALSE)` no longer overcounts when second-order
  bin-averaging is switched on with `second_order_w()`. It was applying the
  prey-bin quadrature twice — once through its bin-averaged prey weight and
  again through the bin-integrated predation kernel that `setPredKernel()`
  builds under `second_order_w` — so the diet was uniformly too large by
  `(1 + beta) / 2`, where `beta` is the grid ratio (9.7% for `NS_params`).
  Summing the diet over prey now reproduces
  `getEncounter() * (1 - getFeedingLevel())` under both schemes, on both the
  FFT and the custom-kernel path. `getDiet(proportion = TRUE)`, the default,
  was unaffected because the factor was uniform and divided out (#474).

- `getTrophicLevel()` had the same quadrature mismatch: its
  trophic-level-weighted numerator was built from the point-sampled
  `getPredKernel()` and a bin-averaged prey weight, while its denominator came
  from `getEncounter()`. Under `second_order_w` the two are no longer the same
  integral, so the reported trophic levels were off by up to 0.06. Numerator
  and denominator now use the same quadrature, and a predator whose prey all
  have trophic level 1 comes out at exactly 2 in both schemes (#474).

- `steady()` now successfully converges when the advective flux scheme is set
  to `"van_leer"` (via `second_order_w`). Previously, the time-stepping
  iteration would fall into a limit cycle because the flux limiter weights
  flipped wildly across cells. We resolved this by introducing an exponential
  moving average relaxation to the limiter `chi` (#522).

- Size-spectrum plots with `size_axis = "l"` now transform number and
  biomass densities from weight to length units using the species-specific
  Jacobian. Biomass density per logarithmic size interval (`power = 2`) uses
  the corresponding logarithmic Jacobian (#469).

- `plotFeedingLevel(include_critical = TRUE)` no longer draws a critical feeding
  level above 1 off the top of the plot. The y axis was fixed to the interval
  from 0 to 1, and now widens when the data need it.

## Documentation

- "Extending mizer" and "Guide: Extending mizer" have been merged into a single
  guide at `guide-extend-mizer`, generated from the `extend-mizer` skill; the old
  address redirects. The topic had been split across three documents that each
  explained `setRateFunction()` in full — the two articles and the "Level 3"
  section of the `change-parameters` skill, which is now a pointer. Each half had
  something the other lacked, and both are now in one place for readers and
  agents alike: the table of the signature and return shape required of every
  replaceable rate, and the rules on respecting the `second_order_w` quadrature
  scheme and on never letting a rate jump as a function of abundance.

  A skill can now keep material out of the agent's copy and in the article, which
  is what lets the worked examples stay evaluated: an `<!-- article-only -->`
  block, the mirror of the existing `<!-- agent-only -->` one. The guide builder
  keeps its content and drops the markers, and `mizerAgents` (>= 0.4.0) drops the
  block as it installs the skill, so a topic still lives in a single file.

- "Using mizer extension packages" is now a guide like the others, generated
  from a new `use-extension-packages` skill, so an agent working in a project
  that loads mizerEcopath, mizerShelf, therMizer or any other extension knows how
  the extension chain works. Being named after its skill, the article moved from
  `using-extension-packages` to `guide-use-extension-packages`; the old address
  redirects. It gains a statement of the two rules that cover almost every
  problem — load the packages before using a model that needs them, and persist
  models with `saveParams()` and `readParams()` rather than `saveRDS()`, since
  the file deliberately holds a base-class object and only `readParams()` puts
  the extension class back.

- The "Point values and bin averages" section of `vignette("numerical_details")`
  now explains where each bin integral is performed and why it must be applied
  exactly once, and a new "The `second_order_w` switch" section documents what
  the flag changes and how to make your own diagnostic second-order accurate.

- New article ["Discontinuous rate functions"](https://sizespectrum.org/mizer/articles/discontinuous_rates.html)
  explains what goes wrong when a custom rate function registered with
  `setRateFunction()` depends discontinuously on the abundances — chattering
  trajectories that keep changing as `dt` is refined, a stalled `steadyNewton()`,
  and an unreliable `getStability()` — why none of the time-stepping methods can
  fix it, and how to avoid it by giving the switch a finite width.

- The topic articles and the AI-agent skills are now one set of documents
  rather than two. Each `inst/skills/<topic>/SKILL.md` is shipped as an agent
  skill and is also the source of the matching `guide-*` article, built by
  `dev_scripts/build_guides.R`. `mizerAgents` (>= 0.3.3) installs the
  skills from the mizer it finds installed, so the two packages can no longer
  drift apart.

  These articles are no longer called cheatsheets. A cheatsheet reminds you of
  something you already know, whereas these assume no prior knowledge, so each
  is now a **guide**. Each is also named after the skill it comes from, and
  titled with that skill's own heading, so a topic has one name instead of
  three: `cheatsheet-fishing` has become `guide-set-up-fishing`, titled "Guide:
  Setting up fishing", from the `set-up-fishing` skill. Every old address
  redirects, but `vignette("cheatsheet-fishing")` does not — see the
  ["Upgrading mizer"](https://sizespectrum.org/mizer/articles/upgrading.html)
  article for the full table of old and new names. The `build-multispecies-model`
  skill was renamed to `build-model` at the same time, since it covers the
  single-species, community and trait constructors too.

  There is now one guide per stage of the workflow, matching the skills one
  for one. Two are new, covering topics that previously had a skill but no
  article: **Running Simulations** (the arguments of `project()`, the four ways
  of giving fishing effort, continuing and comparing runs, and when numerical
  diffusion in the default upwind scheme can damp a real oscillation) and
  **Extending mizer** (a short reference companion to the Extending mizer
  article).

  The former "Model setup and calibration" article has been split into
  **Model setup** and **Steady state and calibration**, which are separate
  tasks reached for at different times. The old URL redirects to the first.

  The remaining guides gain the material that had previously only been
  written on the skill side: the fishing guide now covers `setFishing()`
  and how catchability fixes the units of fishing effort; the calibration
  guide covers `steadyNewton()` and `reproduction_level()`; the model
  setup guide covers saving and reloading a model with
  `saveParams()`/`readParams()`; the changing-parameters guide explains
  that the feeding level is set by `f0` rather than by `h`; and the analysis
  guide recommends `finalParams()` over indexing a time series with
  `idxFinalT()`.

- New `analyse-stability` skill, and the matching "Guide: Dynamic
  Stability" article, cover the experimental stability tools added in this
  version: `getStability()`, `steadyNewton()`,
  `getLimitCycleSim()` and `plotBifurcation()`. Like the other skills it is
  shipped in `inst/skills/` and picked up by `mizerAgents::setup_mizer_agent()`
  from the installed mizer, so an agent's guidance describes the version of
  mizer the project actually runs.

- New `understand-size-spectrum-dynamics` skill, and the matching "Guide:
  Size-Spectrum Dynamics" article, explain how mizer models behave rather than
  which function to call: which quantities you impose and which the model
  produces for itself, the food and predation feedback loops that couple
  species, what sets the slope of the steady-state spectrum and the timescales
  of its dynamics, and the difference between density dependence imposed through
  `R_max` and the density dependence that emerges from the feedback loops. It
  closes with a table mapping a symptom you actually see — a species that
  collapses, oscillates, stops growing before `w_mat`, or refuses to respond to
  fishing — to what to inspect.

- The ["Upgrading mizer"](https://sizespectrum.org/mizer/articles/upgrading.html)
  article is now also shipped as an AI-agent skill, `upgrade-mizer-code`, built
  from `inst/skills/upgrade-mizer-code/SKILL.md` by the same generator as the
  guides. Vignettes are not installed with the package, so an agent helping
  you fix a script that stopped working after an upgrade previously had no
  access to this information and would debug a deliberate, documented change
  from first principles. The skill adds an index that maps the symptom you
  actually see — an "unused argument" error, a deprecation warning, a plot that
  changed, an `identical()` comparison that now fails — to the release that
  caused it and the fix.

# mizer 3.2.1

This patch release fixes how species and gear parameters are handled when they
are assigned. A size given as a weight can now actually be set on a model
specified by lengths, editing one of these parameter tables on its own no
longer validates it on every assignment, and code that sets a species parameter
together with the rate array it determines can now record its change without
triggering a recalculation that would undo it.

## Species parameter changes

- Fixed: `given_species_params<-()` did not apply the length/weight precedence
  rule, so it disagreed with `species_params<-()`, which the documentation
  describes as equivalent apart from its warnings. On a model where a size is
  given both as a weight and as a length, changing the length through
  `given_species_params<-()` left the weight at its old value and warned that
  the length was inconsistent, and it went on warning at every later parameter
  change because the given species parameters were never brought into line.
  The same edit made with `species_params<-()` gave a maturity weight up to 73%
  larger. Both setters now follow the same rule: the one you gave last wins
  (#490). `given_species_params<-()` also preserves the columns of the
  `species_params` slot that mizer does not calculate, which
  `species_params<-()` already did.

- Fixed: `species_params<-()` did not re-derive the calculated species
  parameters that a rate setter owns (`h`, `gamma`, `ks`, `q`, `z0`, `beta`,
  `w_mat25`, …). It rebuilt the species parameters from
  `given_species_params()` but then carried the previously calculated values
  over, so they looked like given values and no longer responded to a change in
  the parameters they are derived from. For example, in a model where `h` is
  derived from the age at maturity, changing `w_mat` left `h` (and with it
  `gamma` and `ks`) at its old value. `species_params<-()` and
  `given_species_params<-()` now agree. Columns that mizer does not calculate
  are still preserved, as are values you supplied yourself.

- `species_params<-()` gains a `recalculate` argument. With
  `species_params(params, recalculate = FALSE) <- value` mizer records the
  values you changed in `given_species_params`, so that they are not calculated
  away later, and stores the parameters you supplied, but does not re-derive
  the calculated species parameters and does not recalculate any of the
  size-dependent rates. This is for code that sets a species parameter together
  with the rate array it determines, where the recalculation would undo the
  caller's own adjustment. Keeping the object consistent is then the caller's
  responsibility.

- New `record_given_species_params()` exports the entry-by-entry change
  detection that `species_params<-()` uses to decide which values to record in
  `given_species_params()`. It is the recording step on its own, without the
  rebuild of the species parameters and the recalculation of all the rates.
  This is for code that has already written into the `species_params` slot
  itself, for example an optimiser that fits a species parameter together with
  the rate array it determines. Such code has to record its changes: a species
  parameter written straight into the slot is silently reverted the next time
  anything triggers a recalculation. If you have a species parameter data frame
  to hand rather than having written into the slot, use
  `species_params(params, recalculate = FALSE) <- value` instead.

- A size that can be given either as a weight or as the length it converts to
  (`w_mat` and `l_mat`, and likewise `w_mat25`, `w_repro_max`, `w_inf`, `w_max`
  and `w_min`) now follows a single rule: **the one given last wins, and if
  both are given at the same time the weight wins**. The other one is set to
  match, so the two never disagree.

  Previously the length always won, whenever and however it had been supplied.
  On a model specified by lengths an assignment like
  `species_params(params)$w_mat[1] <- 100` was therefore silently replaced by
  the value calculated from the unchanged `l_mat`, so weights simply could not
  be set. Setting a length still determines the weight as before.

  This also resolves a contradiction: given both a weight and a length for the
  same size, mizer would report "I will ignore your value for `l_max`" and then
  use `l_max` anyway. It now does what it says.

- The conversion between lengths and weights is now applied per species.
  Previously a single species whose weight and length disagreed caused the
  whole column to be recalculated from the lengths, which also overwrote the
  weights of species whose length was not known at all.

- Mizer now warns, naming the species, when a length and a weight for the same
  size disagree and it changes the length to match the weight. This was
  previously silent and in the other direction: on a model whose stored `w_mat`
  had drifted from `a * l_mat^b`, any assignment to `species_params()` rewrote
  `w_mat` for every species without saying so.

- Editing a `species_params` or `gear_params` data frame on its own no longer
  validates it on every assignment. The `[<-`, `$<-` and `[[<-` methods that
  did this have been removed, so such a data frame now behaves like the plain
  data frame it is: you can leave it in an inconsistent state while you work on
  it, and no conversion, check or warning happens until you validate it, which
  is what `species_params()` does to a plain data frame and what
  `species_params<-()` does to the table it is given. The subclass is still
  preserved, because the base `data.frame` methods preserve it.

  This makes such an assignment about 13 times faster, and it is what lets the
  length/weight rule work as intended: `species_params<-()` can see which
  values you changed relative to the model, so the value you set last wins.
  Note that a data frame edited on its own carries no such history, so if a
  length and a weight both differ they count as given at the same time and the
  weight wins.

- `species_params<-()` no longer errors when the species parameter data frame
  contains a list column (or a column holding S4/other objects), and no longer
  mis-handles a matrix column. The old-vs-new diff now compares such columns
  per species with `identical()` instead of `==`.

- `l2w()` and `w2l()` are about 3 times faster. They were spending about 85% of
  each call on argument checking rather than on the conversion: `assert_that()`
  builds its message whether or not it is needed, and `is(species_params,
  "MizerParams")` consulted the S4 class hierarchy even for a plain data frame.
  They now do the same checks, in the same order and with the same error
  messages, more cheaply. Validating a species parameter data frame is also
  faster, because the length/weight conversion now does the arithmetic inline
  instead of calling these functions once for each of the six size parameters.

- `matchGrowth()` now records the species parameters it scales (`gamma`, `h`,
  `ks` and `k`) among the given species parameters. Previously it wrote them
  straight into the `species_params` slot, so any later species parameter
  change recalculated them from the unscaled given values and silently undid
  part of the match: on `NS_params` the growth rate no longer matched after a
  subsequent assignment to `species_params()`.

## Other changes

- The startup message that `library(mizer)` prints when it first sees a new
  version now only appears when the `major.minor` part of the version changes.
  It links to the release announcement for the series, and there is one of
  those per minor release, so a patch release or a development version no
  longer re-shows an announcement you have already read.

# mizer 3.2.0

This release overhauls how species and resource parameters are set, makes the
extension framework composable regardless of load order, adds a new
`adjustSizeGrid()` function and cheatsheets, and includes a range of smaller
improvements and bug fixes.

For an overview see the 
[release announcement](https://blog.mizer.sizespectrum.org/posts/2026-07-17-mizer-3-2-announcement/)
on the mizer blog.

## Resource setting

- Assigning to `resource_params()` (or one of its components, e.g.
  `resource_params(params)$kappa <- ...`) now immediately rebuilds the
  size-dependent resource rate (`rr_pp`) and capacity (`cc_pp`) arrays from the
  scalars, leaving any manually set (frozen) arrays untouched, exactly as
  `species_params<-` feeds the species rates. As a result changing the
  rate-side scalars `r_pp` or `n` now takes effect (previously the value was
  silently discarded), and successive scalar edits accumulate instead of
  overwriting each other. Assigning to `resource_params()` no longer balances
  the resource; balancing to preserve the steady state is now solely a feature
  of `setResource()`.

- The `resource_rate<-`, `resource_capacity<-`, `resource_level<-` and
  `resource_dynamics<-` setters gained a `balance` argument (default
  unchanged) so it can be switched off, e.g.
  `resource_capacity(params, balance = FALSE) <- my_capacity`.

- `setResource()` no longer silently overwrites a manually set (frozen) rate
  or capacity array when balancing: the frozen array wins and a warning is
  issued. The one exception to the frozen-array protection is that `steady()`
  will rebalance the resource_capacity in order to return a steady state,
  ignoring any freeze.

These changes and how to adapt existing code are described in the new
`vignette("upgrading")` ("Upgrading mizer").

## Species parameter changes

- Modifying species parameters via `species_params<-()` now automatically
  detects your changes, records them in `given_species_params` so they are
  protected from being overwritten by defaults in the future, and silently
  triggers the recalculation of any dependent parameters and rate arrays.
  Previously, `species_params<-()` bypassed the `given_species_params`
  protection and didn't trigger recalculations. This restores expected
  behaviour and makes `species_params<-()` the recommended setter for scripts.

- The `given_species_params<-()` setter remains as an explicit alternative
  that is particularly useful during interactive sessions, because it issues
  warnings if you change a parameter whose effect is overridden by another
  parameter that has already been given.

- Each species parameter default now has a single home: the rate-setting
  function that uses the parameter. `validSpeciesParams()` now only fills in
  defaults for parameters that no single rate setter owns, namely `w_max`,
  `w_repro_max`, `w_mat`, `w_min`, `alpha`, `n`, `a`, `b` and `is_background`.
  The defaults for `p`, `k`, `z_ext`, `d`, `E_ext`, `D_ext` and
  `interaction_resource` are supplied by `setMetabolicRate()`, `setExtMort()`,
  `setExtEncounter()`, `setExtDiffusion()` and `setInteraction()` respectively,
  where they were already being set. Built models are unaffected, because
  `setParams()` calls all the rate-setting functions, but
  `validSpeciesParams()` applied to a bare species parameter data frame now
  returns fewer columns. See the `default_parameters` vignette.

- The `p` argument of `setMetabolicRate()` is deprecated (#459). It never had
  any effect on a `MizerParams` object: such an object always has a `p` column
  already, and the argument was only ever used to fill in a missing one, so it
  was silently ignored. Set the species parameter instead, with
  `species_params(params)$p <- value`. The `p` argument of
  `newMultispeciesParams()` is a different argument and is not affected.

- The default for the metabolic exponent `p` is now `n` rather than `3/4` in
  `setMetabolicRate()`, which is where the default now lives;
  `validSpeciesParams()` no longer sets `p`. No model changes as a result.
  Models built with `newMultispeciesParams()` take `p` from its own `p`
  argument (default `0.7`), which is injected into the species parameter table
  before validation and is untouched by this change, so neither of these
  defaults fires for them. The `validSpeciesParams()` default (`p = n`) only
  ever applied when it was called directly on a bare species parameter data
  frame, which now returns no `p` column, and it shadowed the
  `setMetabolicRate()` default whenever both were in play.

- Default values for the `a` (0.01) and `b` (3) species parameters (for the
  weight-length relationship) are now saved in `species_params` instead of
  being calculated internally by `l2w()` and `w2l()` only when needed.

- The `species_params` data frame is now an S3 subclass of `data.frame`
  (`class = c("species_params", "data.frame")`). It supports class-preserving
  subsetting and subassignment S3 methods, making it safer to use and paving
  the way for future auto-recalculations.

- Columns accessed via `$` on a `species_params` or `gear_params` object now
  return named vectors, where the names are the species names (or "species,
  gear" row names for `gear_params`). For example,
  `species_params(params)$w_mat` now returns a named vector making it easier
  to identify entries. The `species` vector is left unnamed.

- When `sel_func` is set on a `gear_params` object, any argument columns
  required by that selectivity function (other than `w`, `species_params`, and
  `...`) are now automatically added as `NA` columns. This means, for example,
  that setting `gp$sel_func <- "sigmoid_length"` immediately creates the `l25`
  and `l50` columns, ready to be filled in (#431).

- Misspelled column names in the `gear_params` and `species_params` data
  frames are now detected by fuzzy matching against the recognised parameter
  names. A near miss such as `sel_fun` (instead of `sel_func`) triggers a
  warning that suggests the intended name, rather than being silently ignored
  (#442). Columns are only flagged, never renamed, so legitimate custom
  columns are left untouched.

See the new `vignette("upgrading")` ("Upgrading mizer") for how to adapt
existing code to these changes.

## Extension framework

- An installed extension package is now recognised as a dispatching extension
  from the S3 methods it registers for its marker class (e.g.
  `getEncounter.mizerMR`), rather than only from a statically defined S4 marker
  class. This lets extension packages omit the static
  `setClass("mizerFoo", contains = "MizerParams")` and instead let mizer create
  the marker class dynamically when the package is loaded, inserting it at the
  correct place in the S4 hierarchy relative to any other extensions loaded in
  the session. As a result, two independently developed extension packages
  (for example mizerReef and mizerMR) can now be chained in either load order,
  which a static sibling-of-`MizerParams` class prevented.

- `recordExtension()` now prepends a genuinely new extension to the front of the
  object's `@extensions` chain, keeping it ordered outermost-first to match
  `registerExtension()`. Existing entries stay in place.

## New functions

- New `adjustSizeGrid()` function (an S3 generic) adjusts the size grid of
  a `MizerParams` object to a new minimum and/or maximum size. It can both
  expand and truncate (shrink) the grid. For each species it warns if
  truncation discards a non-negligible fraction of the species' biomass, of
  the diet of its smallest individuals, or of the diet of its largest
  individuals.

- Added a `callback` parameter to `project()` to allow user-defined functions
  to be called at each saved time step.

## Other improvements

- Mizer plots no longer produce the unhelpful warning "log-10 transformation
  introduced infinite values" when a logged axis contains zero values (#463).
  
- `setColours()` and `setLinetypes()` now also update the `linecolour` and
  `linetype` entries in `species_params` and `given_species_params` whenever a
  name being set coincides with a species name, so that the choice persists
  with the species rather than only living in the plotting slot.

- `library(mizer)` now prints a one-line startup message the first time you
  load a new mizer version, pointing you to `news(package = "mizer")`. It is
  shown at most once per version and never interrupts a session more than
  that.

- `compareParams()` now checks that the number of size bins, species and gears
  agree before comparing the array-valued slots. When they differ it reports the
  mismatch instead of erroring while trying to compare arrays of incompatible
  dimensions. It also compares the species-parameter tables by matching species
  and parameters by name, so differing species no longer produce a long list of
  per-column length mismatches, and duplicated messages are no longer repeated.

- Error messages that referred to `w_max` as a species' "maximum size" now
  correctly describe it as the upper size-grid boundary, consistent with `w_max`
  being a purely computational parameter.

- `newSingleSpeciesParams()`, `newTraitParams()` and `newCommunityParams()` now
  document why they place the size-grid boundary at the maximum size
  (`w_max = w_repro_max`): because they do not yet set up stochastic growth by
  diffusion, no individual grows beyond `w_repro_max`, so no headroom above it is
  needed. This will be revisited when the constructors gain a diffusion parameter
  (#339).

- `print()` on `ArraySpeciesBySize`, `ArrayTimeBySpecies`, `ArrayResourceBySize`,
  `ArrayTimeByResourceBySize` and `ArrayTimeBySpeciesBySize` objects (as
  returned by `getEncounter()`, `getBiomass()`, `getResourceMort()`, `NResource()`,
  `getFMort()` and similar functions) now prints the array's actual values
  instead of a per-species min/mean/max summary, truncating large arrays to
  fit the console: species are shown as a leading subset, sizes as an evenly
  log-spaced sample across the full size range, and time series as a
  representative sample of time steps that always includes the first and
  last, with a note reporting how much was omitted. A three-dimensional
  `ArrayTimeBySpeciesBySize` object is previewed via its final time slice,
  matching `plot()`'s existing default for that class. Use
  `as.data.frame()` for full, untruncated access to the data.

- `plotYieldGear()` now supports `log_x`, `log_y`, and `log` arguments,
  aligning its arguments with `plotYield()`.

- The upper boundary condition of the size-spectrum solvers now holds the
  abundance at zero above each species' maximum size `w_max`. Without diffusion
  this is automatic and results are unchanged, but with predation diffusion
  switched on it stops a small amount of density leaking to sizes above `w_max`.
  See the "Numerical Details" vignette.

## Documentation

- The "Getting started" guide now includes a self-contained "A worked example:
  the Celtic Sea" section that takes a real ecosystem from raw species parameters
  through building, finding the steady state, calibrating to observed biomasses
  and growth, checking against observed yields, setting the resilience to
  fishing, and projecting a fishing scenario whose sustainable-yield curve is
  interpreted — the whole mizer workflow in one place (#450).

- The first example in the "Getting started" guide no longer prints the
  parameter-default notes from `newMultispeciesParams()`, which were alarming to
  new users out of context (#450).

- Added three new cheatsheets: "Model Setup and Calibration" (building a model,
  finding the steady state, calibrating to observed biomass/yield/growth, and
  projecting), "Fishing" (gears, selectivity functions, catchability, and
  effort), and "Changing Model Parameters" (the distinction between
  `given_species_params()`, `calculated_species_params()` and `species_params()`,
  when changing a species parameter updates a size-dependent rate versus freezing
  it, and how `gear_params()` and the resource setters differ).

- The analysis-and-plotting cheatsheet now covers the newer plotting functions
  (`plotCDF()`, `plotSpectra2()`, `plotSpectraRelative()`, `plotCDF2()`,
  `plot2()`, `plotRelative()`, `animate()`) and the `ArrayResourceBySize`
  class, and corrects the interactive-plot advice for array objects to use
  `plotHover()` (they have no `ggplotly()` method).

- The reference index on the website now opens with an "Overview: the mizer
  workflow" section that frames the whole page as a five-stage pipeline (create →
  calibrate → tune dynamics → project → analyse) with links to the key function
  in each stage, so readers can see how the sections below fit together.

- The reference index on the website now explains the differences between related
  families of functions. It disambiguates the `calibrate...()`/`match...()` and
  `...Biomass`/`...Number` calibration functions and the `plot...()` variants
  (`...2`, `...Relative`, `...ObservedVsModel`), and maps mizer's mortality-rate
  names onto the standard fisheries notation (*M2*, *F*, *Z*).

# mizer 3.1.0

Version 3.1.0 builds on 3.0.0 with an experimental second-order accurate
numerical scheme in size, additional higher-order time-stepping options, and a
range of smaller improvements and bug fixes. Unless you opt in to the
experimental scheme, results are unchanged from 3.0.0.

For an overview see the [release announcement](https://blog.mizer.sizespectrum.org/posts/2026-06-26-mizer-3-1-announcement/) on the mizer blog.

## Experimental second-order in w accuracy

mizer gains an optional, experimental second-order accurate finite-volume scheme
in the size variable `w`. It is controlled by a new `second_order_w` slot and is
switched off by default, so all default results are unchanged (the first-order
path is byte-identical to previous mizer). Enabling it shifts size-integrated
diagnostics and the resource spectrum by `O(Δw)` (more on coarse grids), so
calibrated models may need recalibrating. See `?second_order_w` and the
"Numerical Details" vignette.

- The `second_order_w` slot is a named list with a character entry `flux`
  (`"upwind"`, `"van_leer"` or `"centred"`) selecting the advective
  reconstruction, and a logical entry `bin_average` selecting the bin-averaged
  rate quadratures. A fully second-order scheme needs both. Use the new
  `second_order_w()` / `second_order_w<-()` accessors to get and set them. The
  setter accepts a single logical, which sets both entries — `flux = "van_leer"`
  (a TVD reconstruction that keeps abundances non-negative) and
  `bin_average = TRUE` — or a named vector for individual control, e.g.
  `second_order_w(params) <- c(flux = "centred")` for the unlimited flux that is
  genuinely second order even at extrema. Setting it re-runs `setParams()` to
  rebuild the precomputed arrays. The default `flux = "upwind"`,
  `bin_average = FALSE` is the original first-order upwind scheme. Old objects
  are upgraded automatically.

- `newMultispeciesParams()`, `newTraitParams()`, `newCommunityParams()` and
  `newSingleSpeciesParams()` gain a `second_order_w` argument (default `FALSE`,
  accepting the same values as the setter) that builds the new model with the
  scheme already selected. Under `bin_average` the resource and abundance power
  laws are constructed bin-averaged from the start; the construction-time
  steady-state solve always uses the robust upwind flux, and the chosen flux is
  activated only for the returned model. (#379)

- When `bin_average` is `TRUE`, every point-sampled power law and quadrature that
  feeds the finite-volume update is replaced by its exact bin average over each
  size bin (the bin straddling `w_pp_cutoff` receiving the partial average), so
  the sinks, sources and capacities are consistent with the bin-integrated
  encounter convolution (#374):

  - the external mortality \eqn{z_{ext} w^d} (`setExtMort()`) and external
    diffusion \eqn{D_{ext} w^{n+1}} (`setExtDiffusion()`);
  - the auto-calculated resource intrinsic growth \eqn{r_p w^{n-1}} and carrying
    capacity \eqn{\kappa w^{-\lambda}} (`setResource()`; user-supplied full
    vectors are left untouched), together with the matching initial resource
    abundance \eqn{\kappa w^{-\lambda}} used both for the initial spectrum and
    for the temporary prey spectra behind the default `gamma`/`f0` and consumer
    abundances;
  - the predation kernels, which are now predator- and prey-bin averaged (via
    trapezoid folds), so `getEncounter()`, `getPredRate()`, `getPredMort()`,
    `getResourceMort()` and `getDiffusion()` become second order with no change
    to the rate functions and no extra runtime cost. The predation-diffusion
    integral uses a dedicated Fourier kernel held in the new `ft_pred_kernel_d`
    slot, which carries the extra power of prey size (\eqn{w_p^2 dw_p}, the
    \eqn{\beta^{3s}} Jacobian) that the diffusion integrand needs; in the
    first-order scheme it equals `ft_pred_kernel_e`, so existing models are
    byte-identical. (#384)
  - the gear selectivity (`calc_selectivity()`), so a knife-edge gear gets the
    exact fraction of the straddling bin above the knife edge and fishing
    mortality is second order;
  - the per-capita reproductive investment \eqn{\psi(w) E_r(w)} in `mizerRDI()`
    (the full investment averaged together, not `psi` alone), making
    density-independent reproduction second order.

- The advective growth flux uses the chosen `flux` reconstruction (`"van_leer"`
  or `"centred"`) for a second-order transport step. Combined with the
  bin-averaged diffusion from `getDiffusion()`, the full growth-transport step is
  second order. The diffusion coefficient is consumed from `getDiffusion()`
  directly rather than being re-averaged by the transport routine.

- Size-integrated and size-resolved diagnostics are placed and weighted
  consistently with the finite-volume scheme when `bin_average` is `TRUE`:

  - the summary integrals `getBiomass()`, `getSSB()`, `getYield()`,
    `getYieldGear()`, `getDiet()` and `getTrophicLevel()` use the trapezoidal
    bin-average of the size weight (`getN()` is already exact);
  - the size-resolved bin-average diagnostics — the mortalities `getPredMort()`,
    `getFMort()`, `getMort()`, `getExtMort()` and the reproductive investment
    `getERepro()` — are reported at the geometric bin centre
    \eqn{\sqrt{w_j w_{j+1}}}, the location where a bin average actually lives,
    while point-valued quantities (encounter, growth) stay on the grid nodes.
    The `ArraySpeciesBySize` / `ArrayTimeBySpeciesBySize` classes carry a
    `representation` tag (`"point"` / `"average"`) recording this. (#382)
  - the spectrum plots `plotSpectra()`, `plotlySpectra()`,
    `plotSpectraRelative()` and `animateSpectra()` evaluate both the `w^power`
    weight and the marker location at the geometric bin centre
    \eqn{w^* = w \sqrt{\beta}}, placing each marker as a point on the continuous
    `N w^power` curve instead of misplacing it at the bin edge (the error grows
    with `power`, worst for the common `power = 2` Sheldon plot). (#383)

## Higher-order time-stepping

- `project()` gains a new time-stepping option `method = "tr_bdf2"`. This is an
  L-stable, second-order TR-BDF2 scheme that retains the second-order accuracy of
  `method = "predictor_corrector"` while damping the oscillations the
  Crank-Nicolson corrector can show at large time steps. Like the other methods
  it only requires tridiagonal solves. See the "Numerical Details" vignette.

- Under the second-order methods (`"predictor_corrector"` and `"tr_bdf2"`) the
  resource is now advanced with midpoint resource mortality rather than the
  start-of-step value, and the other components (set via `setComponent()`) now
  also receive a corrector step with the midpoint rates. So the resource and the
  other components are integrated to the same second-order accuracy in time as
  the consumer spectra. The `"euler"` method and the steady states are unchanged.

## Other improvements

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

- `getDiet()` gains a `MizerSim` method and `plotDiet()` for a `MizerSim` now
  accepts a `time_range` argument, computing the diet from the simulated
  abundances at the requested times rather than always using the initial
  abundances. As for the other `MizerSim` plotting functions, `time_range`
  defaults to the final saved time step. When a range spanning several saved
  time steps is given, the consumption rates are averaged over the range and
  then normalised to proportions (rather than averaging the per-step
  proportions, which are normalised independently). (#357)

- New `getFluxGradient()` function returns the flux divergence
  \eqn{(J_{j+1} - J_j)/\Delta w_j} that appears as the transport term in the
  discretised size-spectrum equation. The bin-boundary fluxes are obtained from
  `getFlux()`, so the advective-flux scheme stored in the `flux` entry of the
  `second_order_w` slot is used, with the largest size class closed by the
  boundary condition \eqn{N_{K+1} = 0}. Like `getFlux()`, it has both
  `MizerParams` and `MizerSim` methods, returning an `ArraySpeciesBySize` or
  `ArrayTimeBySpeciesBySize` object respectively.

- `getDiffusion()` now works with a custom predation kernel that depends on
  predator and prey size separately rather than only on their ratio. As for
  `getEncounter()`, when such a kernel has been set (with
  `setPredKernel(params, pred_kernel = ...)`) the diffusion integral is
  evaluated by direct summation over the full predation kernel instead of via
  the FFT method, which assumes a ratio-only kernel. (#373)

- `getTrophicLevel()` and `getTrophicLevelBySpecies()` now assign the resource a
  size-dependent trophic level
  \eqn{T_R(w) = \max(1, 1 + \log(w / w_R) / \log(\beta_R))} instead of treating
  it as trophic level 0. The new `w_R` (average primary-producer size) and
  `beta_R` (average resource predator/prey mass ratio) arguments control this.

- Resource functions now return classed objects that support the same
  convenient `print()`, `summary()`, `plot()` and `as.data.frame()` methods as
  the consumer rate functions. `getResourceMort()`, `initialNResource()`,
  `finalNResource()`, `resource_rate()`, `resource_capacity()` and
  `resource_level()` return an `ArrayResourceBySize` object, and `NResource()`
  returns an `ArrayTimeByResourceBySize` object, so you can now do e.g.
  `plot(getResourceMort(NS_params))` or `plot(NResource(NS_sim))`.

- `summary()` of a `MizerSim` object now reports the fishing effort that was
  actually used during the simulation rather than the model's `initial_effort`.
  Gears whose effort varied over time show the mean effort, flagged with a note
  giving the min-max range.

- Added a new vignette explaining the calculation of default parameter values
  (#189).

## Breaking changes

- The maximum-size species parameters have been clarified (#325). The von
  Bertalanffy asymptotic size `w_inf` is now the required maximum-size parameter
  and is used as the default for `w_repro_max` (previously `w_max`) and `w_mat`.
  `w_max` is now purely a computational boundary (the size grid and plot range)
  and defaults to `1.5 * w_inf`. `w_repro_max` is documented as the size at which
  a typical mature individual invests all its energy into reproduction, not as a
  hard ceiling on size. The default value of the external mortality parameter
  `z0` is now calculated from `w_inf` rather than `w_max`, so the purely
  computational boundary `w_max` no longer influences any model parameter. For
  backwards compatibility, if `w_inf` is not supplied it is taken from
  `w_repro_max` or `w_max`, so existing models and scripts are unaffected;
  however, new models built from the defaults may differ from 3.0.0.

## Bug fixes

- Fixed a bug in `project()` where the abundances of other components (set via
  `setComponent()`) were advanced only once per saved time step instead of once
  per `dt` time step. Their dynamics are now integrated with the same time step
  as the consumer and resource spectra, so results no longer depend on `t_save`.

- `getGrowthCurves()` no longer emits a cryptic
  `collapsing to unique 'x' values` warning when a species' `w_max` lies exactly
  on a grid point. (#413)

- `getRDI()`, `getRDD()` and `getFlux()` on a `MizerSim` object now correctly use
  the simulated time-varying effort instead of the initial effort. (#370)

- `plotCDF()` / `plotlyCDF()` now plot each cumulative value on its bin's
  **upper** edge `w_k + dw_k`, following the inclusive cumulative-sum convention
  (the sum through bin `k` is the integral up to that bin's upper edge). This
  corrects a long-standing one-bin location offset and applies in both the
  default and the second-order schemes; under second-order bin-averaging the CDF
  is then second-order accurate in its placement as well as its increments.
  (#383)

- `plotYield()` now uses `sim2 = NULL` instead of `missing(sim2)` to detect the
  optional second simulation argument, so it works correctly with `do.call()`.
  This is backward-compatible.

- `distanceMaxRelRDI()` now returns `Inf` instead of `NaN` when a previous RDI
  is zero, so `projectToSteady()` no longer mistakes a `NaN` distance for
  convergence.

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
