---
name: upgrade-mizer-code
description: >-
  Diagnose and fix user code or models that broke, changed results, or started
  warning after a mizer upgrade — every documented change from 2.5.4 through 3.3,
  release by release, with the fix. Use whenever a script that "used to work" now
  errors, a deprecation warning appears, plots or numbers differ from a previous
  run, an argument is suddenly unused or rejected (power=, sim=, time_range=,
  setParams() rejecting an argument it does not use), a function has gone
  (matchYields, calibrateYield), a parameter change warns that it cannot take
  effect, an identical() comparison against a saved rate array fails, or `$` on a
  parameter table stopped matching partially. Starts from a symptom index, so
  search it by the message the user actually saw.
---

# Upgrading your mizer code

<!-- agent-only -->

## How to use this skill

Existing model **objects** created with an earlier version are upgraded
automatically when they are loaded, so a saved `MizerParams` or `MizerSim` is
almost never the problem. What changes across releases is *behaviour* and the
*functions the user calls*. Diagnose in this order:

1. **Establish the two versions.** `packageVersion("mizer")` gives the current
   one. Ask the user which version the code last worked with, or infer it from
   the project (a `renv.lock`, a `DESCRIPTION`, the date of the script). This
   fixes the range of releases you have to consider; ignore the rest.
2. **Match the symptom** in the index below rather than debugging from first
   principles. Most upgrade breakages are deliberate, documented changes, and
   reading the model's internals will not reveal that. The index is grouped by
   release: scan only the groups inside the range from step 1.
3. **Read the one section the row names.** The prose for each release lives in
   its own file under `references/`, named in that group's heading. Open only
   the file for the release you matched, and find the `###` heading quoted
   verbatim in the row's Section column. Do not read the other release files.
4. **Apply only the listed fix.** Do not "repair" a model whose numbers moved
   because of a corrected bug — the new numbers are the right ones. Say so, and
   let the user decide whether to recalibrate.

If the symptom is not in the index, it is probably not an upgrade issue at all;
fall back to ordinary debugging, and check `NEWS.md` for the intervening
releases. This skill covers only changes that alter the behaviour of *existing*
code. Purely additive features (new functions, new optional arguments, new
plots) are in the changelog and are not repeated here.

## Symptom index

Each row names a section heading, verbatim, in the reference file for its
release group. `build_guides()` checks that every row resolves to a real
heading and that every heading has at least one row, so a row that does not
match is a bug in this file, not a section you should hunt for elsewhere.

Text in a code span in straight double quotes, `"like this"`, is a literal run
of the message mizer emits, checked against `R/` by the same generator. Search
those first: they are what the user pasted. Messages mizer does not compose
itself — the lifecycle deprecation sentences, base R's `could not find
function` — carry no such quote, so match those rows on the function name.

### mizer 3.3 → 3.4 — `references/mizer-3.4.md`

| Symptom | Cause | Section |
|---|---|---|
| `R_max` or `gamma` no longer reverts to its pre-`scaleModel()` value when the species parameters are recalculated | the rescaled values are now recorded as given | `scaleModel()` records the rescaled parameters as given |
| `given_species_params()` now shows the rescaled `R_max` and `gamma` after `scaleModel()`, `calibrateBiomass()`, `calibrateNumber()`, `matchBiomasses()` or `matchNumbers()` | the rescaling used to bypass the given-value tracking | `scaleModel()` records the rescaled parameters as given |
| Reproduction, or a rate derived from `gamma`, has changed in a model that was rescaled and then had a species parameter edited | the model used to be left with an unscaled `R_max` and `gamma` in an otherwise scaled model | `scaleModel()` records the rescaled parameters as given |

### mizer 3.2 → 3.3 — `references/mizer-3.3.md`

| Symptom | Cause | Section |
|---|---|---|
| `plotSpectra()` or `plotCDF()` errors `"but not contradictory values of both"` | supplying both is no longer silently resolved | `biomass` and `per_log_size` replace `power` |
| A `plotSpectra()` call with both `power` and `biomass`, or any `plotly...()` call with `biomass`, now gives a different plot | `biomass` is no longer ignored | `biomass` and `per_log_size` replace `power` |
| `plotCDF(per_log_size = TRUE)` errors `"A cumulative distribution does not depend on whether the"` | meaningless for a cumulative distribution | `biomass` and `per_log_size` replace `power` |
| `plot(getFluxGradient(...), size_axis = "l")` gives different values, or a `cm^-1/year` label where it used to say `g^-1/year` | the flux gradient is a density and was not recognised as one | Arrays say what kind of value they hold |
| `plot()` of a feeding level, maturity, `psi()` or resource level has a y axis running from 0 to 1 where it used to fit the data | these arrays now declare themselves proportions | Arrays say what kind of value they hold |
| `plot(resource_level(params))` has a linear y axis where it used to be logarithmic | a proportion is plotted on a linear axis by default | Arrays say what kind of value they hold |
| `plotFeedingLevel(include_critical = TRUE)` shows a critical feeding level peak that used to be cut off at 1 | the fixed [0, 1] window is now widened to fit the data | Arrays say what kind of value they hold |
| A custom array plotted with `size_axis = "l"` is not transformed, or is transformed when it should not be | an array declares what it holds with `type` now | Arrays say what kind of value they hold |
| `array_spectrum_power()`, or `spectrum_power =` in an internal plot helper, is no longer found | replaced by the `type` metadata | Arrays say what kind of value they hold |
| New error `"holds a value of type"` from `plot2()` or `plotRelative()` | the two arrays disagree about what they hold, which no pair of axes can carry | Comparison plots transform each array with its own model |
| `plot2(size_axis = "l")` moves the second curve, when the two models differ in `a` or `b` | each array is now transformed with its own weight-length relationship | Comparison plots transform each array with its own model |
| `plotRelative(size_axis = "l")` or `plotSpectraRelative(size_axis = "l")` draws a curve where it used to draw almost nothing | the two series are interpolated onto a common grid instead of matched by equality | Comparison plots transform each array with its own model |
| `plotRelative()` of a density on a length axis gives non-zero differences between models that differ only in `a` or `b` | the two Jacobians no longer cancel out of the ratio | Comparison plots transform each array with its own model |
| `plot2()` or `plotRelative()` now draws a thicker line for `highlight` | the argument used to be swallowed by `...` | Comparison plots transform each array with its own model |
| A `Total` line appears on a plot with `size_axis = "l"` where there used to be none | the total is now summed after the length conversion | The total is summed on the axis it is plotted against |
| `plot(<array>, species = ..., total = TRUE)` gives a bigger total than before | the total is the total of the whole array, not of the selected species | The total is summed on the axis it is plotted against |
| `plotSpectra2()` or `plotSpectraRelative()` with `size_axis = "l"` shows a `Total` line where it used to show none | the total is now formed on the axis being plotted | The total is summed on the axis it is plotted against |
| `plotSpectra2(size_axis = "l", ylim = ..., return_data = TRUE)` returns fewer rows | `ylim` now filters the data there as it does for `plotSpectra()` | The total is summed on the axis it is plotted against |
| `plotBiomass()` or `plot(getBiomass(sim))` draws a background species once where it used to draw it twice | background species are grouped in one pass rather than appended | Background species and non-positive values in time-series plots |
| `plotBiomass(background = FALSE)` no longer shows the background species | it used to skip the appending step and leave them under their own names | Background species and non-positive values in time-series plots |
| `plot(getBiomass(sim), species = ...)` no longer adds the background species | `species` now decides, as it does for the size-spectrum plots | Background species and non-positive values in time-series plots |
| A time-series plot of a model with background species aborts with a base R error about a replacement having 1 row and the data 0 | an empty background group, now handled | Background species and non-positive values in time-series plots |
| `plot(<ArrayTimeBySpecies>, log_y = FALSE)` shows zero or negative values it used to drop | the `1e-20` floor belongs to a logarithmic axis only | Background species and non-positive values in time-series plots |
| `animate(NResource(sim), size_axis = "l")` labels its y axis `1/cm` where it used to say `1/g` | the label now follows the size axis, as the values already did | Animations follow the array type and the size axis |
| `animate(getFeedingLevel(sim))` has a linear y axis from 0 to 1 where it used to be logarithmic and fitted to the data | an animation of a proportion is scaled like a static plot of one | Animations follow the array type and the size axis |
| `summary()` of a rate array gives a much smaller `Max` (or larger `Min`) than before | it now covers each species' own size range, as `plot()` always has | `summary()` of an array covers the same sizes as `plot()` |
| `summary()` of an array reports `NA` where it used to report `Inf` or `-Inf` | an empty selection is now reported as missing rather than reduced | `summary()` of an array covers the same sizes as `plot()` |
| `getProportionOfLargeFish(params)` gives a different value, or no longer differs from the `MizerSim` value | weights were recycled down the columns for all but the first species | `getProportionOfLargeFish()` on a `MizerParams` object was wrong |
| `plotYieldObservedVsModel()` errors `"You have not provided values for the column 'yield_observed'"`, but `gear_params()` has the column | the plot used to read only the species parameters | `yield_observed` belongs to the gear parameters |
| Setting a weight on a length-based model no longer gets undone | length/weight precedence: the one given last wins | Length and weight parameters follow the one you gave last |
| `given_species_params<-()` now changes a weight when you change its length | both setters apply the same precedence rule | Length and weight parameters follow the one you gave last |
| Repeated `"is not consistent with the value of"` warning has stopped | the given species parameters are brought into line | Length and weight parameters follow the one you gave last |
| Warning `"marking it as missing so that its default will"` be used, where the message used to say the value was set to `NA` | the message now names the value that is actually stored | An invalid `w_mat25` is replaced by its default |
| Maturity, or a rate derived from it, has changed in a model whose species parameters carry `l_mat25` | a rejected `w_mat25` is no longer restored from its length, so the default is used | An invalid `w_mat25` is replaced by its default |
| New message `"is set by length, but"` from `setFishing()`, `calc_selectivity()` or a model build | a length-based gear is converting lengths with defaulted weight-length parameters | Mizer says when it defaults the weight-length parameters |
| New message `"column so using a = 0.01 in w = a l^b"` or `"the isometric default b = 3"` when a model is built | mizer now reports the weight-length defaults it fills in | Mizer says when it defaults the weight-length parameters |
| New warning `"has not taken effect because the"` after a species or resource parameter change | the rate it feeds was set by hand and is no longer calculated | Species parameter setters distinguish edits from declarations |
| `given_species_params<-()` warning behaviour changed: warnings appear only there, clearing a given value can warn, and adding an all-`NA` column does not | it reports actual changes to the authoritative given-value table | Species parameter setters distinguish edits from declarations |
| Marking a current calculated value as given, or changing an observation or custom column, no longer rebuilds all rate arrays | provenance-only and uncached changes need no rebuild | Species parameter setters distinguish edits from declarations |
| `given_species_params()` loses defaults such as `a` and `b`, or a derived value starts moving again | validation-filled defaults are no longer mistaken for user input | Species parameter setters distinguish edits from declarations |
| Message `"I have removed the species parameter column"` | a column missing from an assigned species parameter table is now withdrawn | A column dropped from an assigned species parameter table is removed |
| Species parameters the code never touched have moved, or been dropped, after assigning a table built from a few columns | the columns left out of that table counted as withdrawn | A column dropped from an assigned species parameter table is removed |
| `given_species_params(params)$gamma <- NULL` now changes `species_params(params)$gamma` | a removal hands the parameter back to mizer's calculation and rebuilds | A column dropped from an assigned species parameter table is removed |
| A custom species parameter column disappears where it used to survive, or `calculated_species_params()` no longer reports it | mizer cannot recalculate a column of your own, so withdrawing it removes it | A column dropped from an assigned species parameter table is removed |
| A column read with `$` is now `NULL`, warning `"Earlier versions of mizer partially matched the column name"` | `$` no longer partially matches column names | `$` on a parameter table no longer partially matches |
| Length-weight conversion via `$a`/`$b` gave `alpha`/`beta` values | partial matching, now fixed | `$` on a parameter table no longer partially matches |
| Size grid or results changed in a model built with a small `min_w` | `w_min` is no longer reset to 0.001 | `w_min` survives a rebuild of the species parameters |
| A recalculated `gamma` or `f0` is wildly different, in a model whose `search_vol` was set by hand | the frozen array used to block mizer's own unit-gamma calculation | The defaults for `gamma` and `f0` are measured on mizer's own reference state |
| On a model using an extension that changes the encounter rate, `species_params(params)$gamma` moves by the same factor every time the species parameters are touched | the `gamma` default used to be measured through the extension's `projectEncounter()` method | The defaults for `gamma` and `f0` are measured on mizer's own reference state |
| `gamma` or `f0` changes on a model with an encounter function registered by `setRateFunction()` | that function no longer enters the calculation of these two defaults | The defaults for `gamma` and `f0` are measured on mizer's own reference state |
| A recalculated `f0` is lower on a model with external encounter, `other_encounter()` or a component `encounter_fun`, or `f0` → `gamma` → `f0` now round-trips | additive encounter contributions are no longer part of the power-law reference state | The defaults for `gamma` and `f0` are measured on mizer's own reference state |
| Error ``"Could not calculate a default `gamma` for the following species:"`` from `species_params<-()` or `upgradeParams()` on an extension model | the extension zeroed the encounter rate for a species while the `gamma` default was being measured through it | The defaults for `gamma` and `f0` are measured on mizer's own reference state |
| The same error on a model that builds a species with `interaction_resource = 0` and an external, free-standing or component encounter | that additive contribution used to stand in for the missing predation encounter and yielded a `gamma` many orders of magnitude too large | The defaults for `gamma` and `f0` are measured on mizer's own reference state |
| Setting `f0 = 1` errors `"must be finite and in the interval [0, 1)"`, even though `gamma` is supplied | every supplied target feeding level is now validated | `f0` is always validated |
| `gamma`, `q` or feeding levels change after setting resource `kappa` or `lambda` | calculated search-volume parameters now follow the resource power law | Resource scalars refresh calculated `gamma` and `q` |
| `setExtMort(z0pre = ...)`, `setExtMort(z0exp = ...)` or the same arguments to `setParams()` now warn ``"`z0` is already present in `given_species_params` for every species"`` | `z0` was already present in `given_species_params()` for every species, so the arguments were ignored | `setExtMort()` warns when `z0pre` or `z0exp` is ignored |
| A message that used to appear no longer does, with `info_level = 0` | `info_level = 0` now silences everything | One report, one switch |
| `expect_message()` on a mizer call fails, or `options(warn = 2)` trips | reports are warnings where they were messages, and are collected until the end of the call | One report, one switch |
| Warning `"very close to standard parameter names"` appears once where it used to appear many times, or a warning count in a test has changed | the misspelling check now runs once, where a column enters the model | One report, one switch |
| A misspelled species parameter column is no longer flagged at all | the report now follows `info_level`, and the model was built with `info_level = 0` | One report, one switch |
| R errors with `could not find function`, naming `steadyNewton` | it never shipped; the Newton solver is now the `solver = "newton"` argument of the two finders | The steady-state finders have new names |
| A code review, a linter or the reference index calls `steady()` or `projectToSteady()` superseded | both were renamed for what they keep fixed | The steady-state finders have new names |
| `projectToSteady(return_sim = TRUE)` has no equivalent on the new functions | the return type no longer depends on an argument | The steady-state finders have new names |
| `conv$type` or `conv$settled` on a `"convergence"` attribute is `NULL`, or an `expect_named()` on it fails | the attribute now has `termination`, `converged` and `attractor` | The convergence attribute has a new shape |
| A limit cycle used to be reported as a converged steady state | the fixed-point claim is now made on the measured biomass drift | The convergence attribute has a new shape |
| `attr(params, "convergence")$termination` is `"distance_tolerance"` after `steady()` | the run stopped on the distance criterion with the model still drifting | The convergence attribute has a new shape |
| Message `"which is below the distance tolerance, but the"`, or a `tuneSteadyState()` run that goes to `t_max` where `steady()` stopped early | the new finders also require the biomasses to have stopped moving | The convergence attribute has a new shape |
| `steady()` says `"Reached the convergence tolerance after"` where it used to announce that convergence was achieved | the message says what was actually tested | `steady()` reports the tolerance it reached rather than announcing convergence |
| An `expect_message()` or a grep for the wording of the `steady()` success message no longer matches | the message was reworded | `steady()` reports the tolerance it reached rather than announcing convergence |
| `steady()` adds `"Reduce the tolerance on the distance function to converge further."` to its convergence message | the state reached is not a fixed point | `steady()` reports the tolerance it reached rather than announcing convergence |
| `steady()` gives different results for a model with seasonal or otherwise time-dependent rates | every block used to restart the clock at zero | The steady-state run advances time like `project()` does |
| `projectToSteady()` finds a limit cycle much earlier, or one it used to miss | ignores the first half of the simulation | `projectToSteady()` ignores initial transients |
| `projectUntilSettled()` or `tuneSteadyState()` converges where it used to run to `t_max`, or `distanceSSLogN()` returns a smaller number | size classes holding a negligible share of a species' biomass are no longer counted | A size class holding no fish no longer blocks convergence |
| Message `"reached is a fixed point: the biomasses change at only"` | a run stopped at `t_max` on a state that is a fixed point all the same | A size class holding no fish no longer blocks convergence |
| New warning `"and stability machinery covers the consumers and the resource only"` | a component with its own dynamics is held fixed by these tools | The steady-state tools hold other components fixed |
| New warning `"not be rebalanced and the preserved resource abundance need not be a steady"` from `steady()` | a custom `resource_dynamics` with no `balance_<dynamics>()` function | The steady-state tools hold other components fixed |
| `summary(params)` has an extra `"Steady state:"` block | new steadiness verdict | `summary()` reports the steady state |
| A `match…()` call now reports `"has rescaled the model and so moved it off its steady state"` | the functions say so themselves now | The `match…()` functions announce that they broke the steady state |
| `matchNumbers()` no longer says that it rescaled the model, or no longer updates `time_modified` | it now returns early when it has nothing to match | The `match…()` functions announce that they broke the steady state |
| `getReproductionLevel()` changed after a `matchNumbers()` call that matched nothing | that call no longer re-tunes the reproduction parameters | The `match…()` functions announce that they broke the steady state |
| Absolute diet values dropped ~10%, only with `second_order_w` bin-averaging | double-counted prey-bin quadrature, now fixed | Fixes under the second-order size scheme |
| Trophic levels moved slightly, only with `second_order_w` | numerator and denominator now use one quadrature | Fixes under the second-order size scheme |
| `getN()` over a size range moved slightly, only with `second_order_w` | the size-range window is now bin-averaged too | Fixes under the second-order size scheme |
| A model calibrated or matched to observations moved slightly, only with `second_order_w` | the calibration functions hand-rolled a first-order sum and now use the model's own quadrature | Fixes under the second-order size scheme |
| `plotBiomassObservedVsModel()` shows a matched species off the 1:1 line, only with `second_order_w` | the plot hand-rolled its own biomass integral | Fixes under the second-order size scheme |
| `plotYieldObservedVsModel()` model yields rose by 10-20%, or a model looked like it under-predicted yields while `getYield()` disagreed | the plot hand-rolled a first-order yield integral instead of calling `getYield()` | Fixes under the second-order size scheme |
| `plot(getFMort(sim))` draws at other sizes than before, only with `second_order_w` bin-averaging | the time slice kept its `representation`, so bin averages sit at the bin centres | Fixes under the second-order size scheme |
| `steady()` now converges on a `van_leer` model where it used to report a limit cycle, or never settle | the flux limiter is relaxed between iterations | Fixes under the second-order size scheme |
| `setParams()` or `setResource()` errors ``"`setParams()` does not have"`` | unknown arguments are no longer silently ignored | `setParams()` rejects arguments it does not use |
| A `setParams(resource_rate = )` / `setParams(kappa = )` call that ran fine now errors | `setParams()` never set the resource; use `setResource()` | `setParams()` rejects arguments it does not use |
| A resource change made via `setParams()` never showed up in the model | the argument was silently dropped in every version before 3.3 | `setParams()` rejects arguments it does not use |
| Deprecation warning for `r_pp`/`kappa` now names `setResource()` instead of `setParams()` | the old recommendation pointed at a no-op | `setParams()` rejects arguments it does not use |
| Old code calls `getCatchability()`, `getPredKernel()`, `getMetabolicRate()` or another `get`-prefixed array accessor, and the user asks whether it still works | it does, silently and permanently; the bare name is the one to use in new code | One name for each stored rate array |
| `getInteraction()`, `getResourceRate()`, `getResourceLevel()`, `getResourceCapacity()` or `getResourceDynamics()` used to warn and now does not | each became a plain alias of its bare name along with the other `get`-prefixed accessors | One name for each stored rate array |
| `getExtMort.MizerParams`, `getInteraction.MizerParams` and friends no longer found as S3 methods | the `get` names are aliases of the bare names now; dispatch happens on the bare name | One name for each stored rate array |
| R errors with `could not find function`, naming `matchYields` or `calibrateYield` | both removed after deprecation in 2.6.0 | `matchYields()` and `calibrateYield()` have been removed |
| `compareParams()` now reports differences it used to miss | relative tolerance for species parameters | `compareParams()` compares small parameters properly |
| After `rm(list = ls())`, or in the second and later examples of an extension package's `R CMD check`, mizer reports that the extension's class is not a defined class | the dynamic marker classes lived in `.GlobalEnv` and were wiped along with it | Extension marker classes are created and repaired by mizer |
| After `devtools::load_all()` on an extension package, `coerceToExtensionClass()` reports no method or default for coercing `MizerParams` to the extension class | repeated registration left the chain's dynamic marker classes unrepaired | Extension marker classes are created and repaired by mizer |
| `setComponent()` errors `"already a rate contribution registered under the name"` | the name is already taken by a free-standing `other_mort()` or `other_encounter()` entry, which the component used to take over silently | `other_mort()` and `other_encounter()` register contributions that have no component |
| A component's `encounter_fun` without `...` errors about an unused `t` argument | encounter contributions now receive the current simulation time, matching mortality contributions | `other_mort()` and `other_encounter()` register contributions that have no component |
| `vignette("cheatsheet-fishing")` (or any other `cheatsheet-…`) finds nothing | the cheatsheet articles were renamed after the skills they come from | The cheatsheet articles are now called guides |
| The `build-multispecies-model` skill is not found | renamed to **build-model** | The cheatsheet articles are now called guides |
| A link to `articles/using-extension-packages.html` | renamed to **guide-use-extension-packages** when it became a generated guide; the old address redirects | The cheatsheet articles are now called guides |
| A link to `articles/extending-mizer.html` | merged into **guide-extend-mizer**, which was previously a separate shorter guide; the old address redirects | The cheatsheet articles are now called guides |
| A link to `articles/creating-extension-packages.html` | renamed to **guide-create-extension-package** when it became a generated guide; the old address redirects | The cheatsheet articles are now called guides |

### mizer 3.1 → 3.2 — `references/mizer-3.2.md`

| Symptom | Cause | Section |
|---|---|---|
| Custom column written with `species_params<-()` now survives, or `w_mat`/`w_max` recalculate unbidden | `species_params<-()` diffs and protects | `species_params<-()` now detects and protects changes |
| `cc_pp` / `rr_pp` change as soon as a resource scalar is set | resource assignment rebuilds the arrays | Setting resource parameters |
| Code that set a resource scalar and then called `setResource()` to apply it | the assignment already rebuilt `cc_pp`/`rr_pp`, so the call is redundant | Assigning to `resource_params()` now updates the resource arrays |
| Resource steady state shifts after setting `kappa` or `r_pp` | assignment does not balance | Assigning to `resource_params()` does not balance the resource |
| Warning `"has been set manually and so it was"` not rebalanced | frozen arrays are protected | Frozen arrays are protected from incidental balancing |
| Setting one resource array silently rebalanced the other, and you wanted it left alone | balancing is now optional | The resource setters gained a `balance` argument |
| An extension's `setClass("mizerFoo", contains = "MizerParams")` is no longer needed, or two extensions now chain in either load order | marker classes are created dynamically from the registered S3 methods | Extension packages: dynamic marker classes |
| `class(species_params(params))` no longer returns just `data.frame` | it is an S3 subclass now | The `species_params` data frame is now an S3 subclass |
| A column extracted with `$` is unexpectedly named | named by species | Accessing a column with `$` now returns a named vector |
| `gear_params` has extra `NA` columns after setting `sel_func` | argument columns added automatically | Setting `sel_func` adds the required argument columns |
| `species_params(df)` on a bare data frame adds columns or errors | full validation now runs | Passing a data frame to `species_params()` / `given_species_params()` now validates it |
| Printing a rate getter's result is truncated | new `print()` methods | Printing of mizer array objects shows the values |
| Density above `w_max` disappeared, diffusion switched on | new upper boundary condition | Upper boundary condition at `w_max` |

### mizer 3.0 → 3.1 — `references/mizer-3.1.md`

| Symptom | Cause | Section |
|---|---|---|
| New model built from defaults differs from 3.0 | `w_inf` is now the primary maximum size | Maximum-size species parameters clarified |
| Trophic levels are higher than before | resource has a size-dependent trophic level | `getTrophicLevel()` gives the resource a size-dependent trophic level |
| `summary()` of a `MizerSim` reports a different effort | reports effort actually used | Bug fixes that change results |
| Results change with `t_save` in a model with other components | other components now advanced every `dt` | Bug fixes that change results |
| `plotCDF()` curves shifted by one bin | bin placement corrected | Bug fixes that change results |
| `projectToSteady()` converges differently | `distanceMaxRelRDI()` returns `Inf`, not `NaN` | Bug fixes that change results |
| `method = "predictor_corrector"` results differ slightly | resource advanced at the midpoint | Second-order methods advance the resource at the midpoint |
| Size-integrated diagnostics shifted after enabling `second_order_w` | opt-in scheme change | Opting in to the second-order-in-size scheme |

### mizer 2.5.4 → 3.0 — `references/mizer-3.0.md`

| Symptom | Cause | Section |
|---|---|---|
| `unused argument (sim = ...)` from `plotBiomass()`, `plotYield()`, `plotYieldGear()` | first argument renamed to `object` | Renamed arguments and changed defaults (breaking changes) |
| `unused argument (time_range = ...)` from `plotDiet()` | removed in 3.0, back for `MizerSim` in 3.1 | Renamed arguments and changed defaults (breaking changes) |
| `setInitialValues()` warns that it is deprecated | replaced by `finalParams()` | `setInitialValues()` is deprecated |
| `identical()` against a stored rate array is now `FALSE`, values unchanged | rate getters return classed arrays | Rate getters return classed array objects |
| Code indexing `getMort()` / `getPredRate()` by dimname fails | dimnames are now `sp` and `w` | Renamed arguments and changed defaults (breaking changes) |
| `selectMethod()` / `getMethod()` on `plot` or `summary` fails | they are S3 methods now | `plot()` and `summary()` are now S3 methods |
| Growth of large or food-limited individuals no longer negative | growth clamped at zero | Growth can no longer be negative |
| Growth, mortality or encounter changed after setting `use_predation_diffusion`, `z_ext`, `d`, `E_ext` or `D_ext` | new opt-in terms, all defaulting to leave the model unchanged | Predation diffusion is available but off by default |
| `project()` warns `"Appending a simulation run with dt ="` or `"Appending a simulation run with method ="` | inherited from `sim_params` | `project()` timing and effort handling |
| Simulation has one more saved time step than before | state at `t_max` is always saved | `project()` timing and effort handling |
| Simulation length or save times changed when an effort array is given | `t_max`/`t_save` now respected | `project()` timing and effort handling |
| `plotBiomassObservedVsModel()` no longer shows ratios | default is now `ratio = FALSE` | Renamed arguments and changed defaults (breaking changes) |
| `plotSpectra()` axes look different | y-limit auto-scales, new default lower size limit | Bug fixes that change results |
For guidance on which accessor to reach for once the diagnosis is made, see the
`change-parameters` skill.

<!-- /agent-only -->

This article collects the changes that may require you to update your own code
or models when you upgrade to a new release of mizer. Existing model objects
created with an earlier version continue to load and run — they are upgraded
automatically — so the notes below are about changes in *behaviour* and in the
functions you call, not about stored objects. The changes are grouped by the
release in which they took effect, most recent first.

Only changes that can alter the behaviour of *existing* code are listed. The
many purely additive features (new functions, new optional arguments, new
plots) are described in the [changelog](https://sizespectrum.org/mizer/news/index.html)
and are not repeated here.
