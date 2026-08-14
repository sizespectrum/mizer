---
name: upgrade-mizer-code
description: >-
  Diagnose and fix user code or models that broke, changed results, or started
  warning after a mizer upgrade. Use whenever a script that "used to work" now
  errors, a deprecation warning appears, plots or numbers differ from a previous
  run, an argument is suddenly "unused" (sim=, time_range=), an identical()
  comparison against a saved rate array fails, or the user asks how to move code
  from mizer 2.5/3.0/3.1 to a later version. Lists every documented change in
  behaviour, release by release, with the fix.
---

# Upgrading mizer

<!-- agent-only -->

## How to use this skill

Existing model **objects** created with an earlier version are upgraded
automatically when they are loaded, so a saved `MizerParams` or `MizerSim` is
almost never the problem. What changes across releases is *behaviour* and the
*functions the user calls*. Diagnose in this order:

1. **Establish the two versions.** `packageVersion("mizer")` gives the current
   one. Ask the user which version the code last worked with, or infer it from
   the project (a `renv.lock`, a `DESCRIPTION`, the date of the script). Every
   section below is keyed by release.
2. **Match the symptom** in the table below rather than debugging from first
   principles. Most upgrade breakages are deliberate, documented changes, and
   reading the model's internals will not reveal that.
3. **Apply only the listed fix.** Do not "repair" a model whose numbers moved
   because of a corrected bug — the new numbers are the right ones. Say so, and
   let the user decide whether to recalibrate.

If the symptom is not in the table, it is probably not an upgrade issue at all;
fall back to ordinary debugging, and check `NEWS.md` for the intervening
releases. This skill covers only changes that alter the behaviour of *existing*
code. Purely additive features (new functions, new optional arguments, new
plots) are in the changelog and are not repeated here.

## Symptom index

| Symptom | Cause | Section |
|---|---|---|
| `plotSpectra()` or `plotCDF()` errors that `power` and `biomass` are contradictory | supplying both is no longer silently resolved | `biomass` and `per_log_size` replace `power` (3.3) |
| A `plotSpectra()` call with both `power` and `biomass`, or any `plotly...()` call with `biomass`, now gives a different plot | `biomass` is no longer ignored | `biomass` and `per_log_size` replace `power` (3.3) |
| `plotCDF(per_log_size = TRUE)` errors | meaningless for a cumulative distribution | `biomass` and `per_log_size` replace `power` (3.3) |
| New warning that a change to a species or resource parameter "has not taken effect" | the rate it feeds was set by hand and is no longer calculated | A change that cannot take effect now warns (3.3) |
| A message that used to appear no longer does, with `info_level = 0` | `info_level = 0` now silences everything | One report, one switch (3.3) |
| `expect_message()` on a mizer call fails, or `options(warn = 2)` trips | reports are warnings where they were messages, and are collected until the end of the call | One report, one switch (3.3) |
| A column read with `$` is now `NULL`, with a warning naming another column | `$` no longer partially matches column names | `$` no longer partially matches (3.3) |
| Length-weight conversion via `$a`/`$b` gave `alpha`/`beta` values | partial matching, now fixed | `$` no longer partially matches (3.3) |
| Absolute diet values dropped ~10%, only with `second_order_w` bin-averaging | double-counted prey-bin quadrature, now fixed | Quadrature fixes under `second_order_w` (3.3) |
| Trophic levels moved slightly, only with `second_order_w` | numerator and denominator now use one quadrature | Quadrature fixes under `second_order_w` (3.3) |
| `getProportionOfLargeFish(params)` gives a different value, or no longer differs from the `MizerSim` value | weights were recycled down the columns for all but the first species | `getProportionOfLargeFish()` on a `MizerParams` object was wrong (3.3) |
| `getN()` over a size range moved slightly, only with `second_order_w` | the size-range window is now bin-averaged too | `getN()` respects the quadrature scheme at the ends of the size range (3.3) |
| Setting a weight on a length-based model no longer gets undone | length/weight precedence: the one given last wins | Length and weight follow the one you gave last (3.3) |
| `given_species_params<-()` now changes a weight when you change its length | both setters apply the same precedence rule | Length and weight follow the one you gave last (3.3) |
| Repeated "`l_mat` is not consistent with `w_mat`" warning has stopped | the given species parameters are brought into line | Length and weight follow the one you gave last (3.3) |
| Size grid or results changed in a model built with a small `min_w` | `w_min` is no longer reset to 0.001 | `w_min` survives a rebuild (3.3) |
| `compareParams()` now reports differences it used to miss | relative tolerance for species parameters | `compareParams()` compares small parameters (3.3) |
| A recalculated `gamma` or `f0` is wildly different, in a model whose `search_vol` was set by hand | the frozen array used to block mizer's own unit-gamma calculation | Defaults for `gamma` and `f0` ignore a hand-set search volume (3.3) |
| `setParams()` or `setResource()` errors that it "does not have an argument" | unknown arguments are no longer silently ignored | `setParams()` rejects arguments it does not use (3.3) |
| A `setParams(resource_rate = )` / `setParams(kappa = )` call that ran fine now errors | `setParams()` never set the resource; use `setResource()` | `setParams()` rejects arguments it does not use (3.3) |
| A resource change made via `setParams()` never showed up in the model | the argument was silently dropped in every version before 3.3 | `setParams()` rejects arguments it does not use (3.3) |
| Deprecation warning for `r_pp`/`kappa` now names `setResource()` instead of `setParams()` | the old recommendation pointed at a no-op | `setParams()` rejects arguments it does not use (3.3) |
| Deprecation warning from `getCatchability()`, `getPredKernel()`, `getMetabolicRate()` or another `get`-prefixed array accessor | the bare name is now the only supported one | One name for each stored rate array (3.3) |
| `getExtMort.MizerParams` and friends no longer found as S3 methods | the `get` names are plain forwarding functions now; dispatch happens on the bare name | One name for each stored rate array (3.3) |
| `unused argument (sim = ...)` from `plotBiomass()`, `plotYield()`, `plotYieldGear()` | first argument renamed to `object` | Renamed arguments and changed defaults (3.0) |
| `unused argument (time_range = ...)` from `plotDiet()` | removed in 3.0, back for `MizerSim` in 3.1 | Renamed arguments and changed defaults (3.0) |
| `setInitialValues()` warns that it is deprecated | replaced by `finalParams()` | `setInitialValues()` is deprecated (3.0) |
| `identical()` against a stored rate array is now `FALSE`, values unchanged | rate getters return classed arrays | Rate getters return classed array objects (3.0) |
| Code indexing `getMort()` / `getPredRate()` by dimname fails | dimnames are now `sp` and `w` | Renamed arguments and changed defaults (3.0) |
| `selectMethod()` / `getMethod()` on `plot` or `summary` fails | they are S3 methods now | `plot()` and `summary()` are now S3 methods (3.0) |
| Growth of large or food-limited individuals no longer negative | growth clamped at zero | Growth can no longer be negative (3.0) |
| `project()` warns that `dt` or `method` differ from the stored ones | inherited from `sim_params` | `project()` timing and effort handling (3.0) |
| Simulation has one more saved time step than before | state at `t_max` is always saved | `project()` timing and effort handling (3.0) |
| Simulation length or save times changed when an effort array is given | `t_max`/`t_save` now respected | `project()` timing and effort handling (3.0) |
| `plotBiomassObservedVsModel()` no longer shows ratios | default is now `ratio = FALSE` | Renamed arguments and changed defaults (3.0) |
| `plotSpectra()` axes look different | y-limit auto-scales, new default lower size limit | Bug fixes that change results (3.0) |
| New model built from defaults differs from 3.0 | `w_inf` is now the primary maximum size | Maximum-size species parameters clarified (3.1) |
| Trophic levels are higher than before | resource has a size-dependent trophic level | `getTrophicLevel()` (3.1) |
| `summary()` of a `MizerSim` reports a different effort | reports effort actually used | Bug fixes that change results (3.1) |
| Results change with `t_save` in a model with other components | other components now advanced every `dt` | Bug fixes that change results (3.1) |
| `plotCDF()` curves shifted by one bin | bin placement corrected | Bug fixes that change results (3.1) |
| `projectToSteady()` converges differently | `distanceMaxRelRDI()` returns `Inf`, not `NaN` | Bug fixes that change results (3.1) |
| `method = "predictor_corrector"` results differ slightly | resource advanced at the midpoint | Second-order methods (3.1) |
| Size-integrated diagnostics shifted after enabling `second_order_w` | opt-in scheme change | Opting in to the second-order-in-size scheme (3.1) |
| Custom column written with `species_params<-()` now survives, or `w_mat`/`w_max` recalculate unbidden | `species_params<-()` diffs and protects | `species_params<-()` detects and protects changes (3.2) |
| `cc_pp` / `rr_pp` change as soon as a resource scalar is set | resource assignment rebuilds the arrays | Setting resource parameters (3.2) |
| Resource steady state shifts after setting `kappa` or `r_pp` | assignment does not balance | Assigning to `resource_params()` does not balance (3.2) |
| Warning that a manually set resource array was kept | frozen arrays are protected | Frozen arrays are protected (3.2) |
| `class(species_params(params))` is not `"data.frame"` | it is an S3 subclass now | The `species_params` data frame is an S3 subclass (3.2) |
| A column extracted with `$` is unexpectedly named | named by species | Accessing a column with `$` returns a named vector (3.2) |
| `gear_params` has extra `NA` columns after setting `sel_func` | argument columns added automatically | Setting `sel_func` adds the required argument columns (3.2) |
| `species_params(df)` on a bare data frame adds columns or errors | full validation now runs | Passing a data frame now validates it (3.2) |
| Printing a rate getter's result is truncated | new `print()` methods | Printing of mizer array objects (3.2) |
| Density above `w_max` disappeared, diffusion switched on | new upper boundary condition | Upper boundary condition at `w_max` (3.2) |

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

## Upgrading from mizer 3.2 to 3.3

Most of the changes in this release are corrections. Results move only for
models that had opted in to second-order bin-averaging, that set `min_w` below
the default, or that specify sizes as lengths. The one change to an interface is
in the spectrum plots.

### `biomass` and `per_log_size` replace `power`

`plotSpectra()`, `plotSpectra2()`, `plotCDF()`, `plotCDF2()` and `animate()`
now describe the plotted quantity with two independent arguments: `biomass`
chooses a biomass density rather than a number density, and the new
`per_log_size` chooses a density with respect to logarithmic size rather than
with respect to size. The `power` of the weight multiplying the number density
is the sum of the two:

| | `per_log_size = FALSE` | `per_log_size = TRUE` |
|---|---|---|
| `biomass = FALSE` | `power = 0` | `power = 1` |
| `biomass = TRUE`  | `power = 1` | `power = 2` |

`power` still works and is still the only way to ask for a power that is not
the sum of the two flags, so calls that pass only `power` are unaffected. Two
things change:

- Passing `power` together with `biomass` used to ignore `biomass` silently.
  Now the two must agree, or you get an error (#501). Where they do agree the
  call is honoured: `plotSpectra(sim, power = 1, biomass = FALSE)` used to plot
  the biomass density, and now plots the number density with respect to
  logarithmic size — the same numbers, but labelled correctly and, with
  `size_axis = "l"`, converted to a length axis with the logarithmic Jacobian
  rather than the density one. If you meant the biomass density, drop the
  `biomass` argument. The same applies to `plotlySpectra()`, `plotlyCDF()`,
  `plotlySpectra2()` and `plotlyCDF2()`, which passed `power` on internally and
  so ignored `biomass` even when you gave only `biomass`: those calls now plot
  what they were asked for.
- `plotCDF()` and `plotCDF2()` do not accept `per_log_size`, because
  integrating a density over size gives the same cumulative quantity either
  way. Use `biomass` on its own there.

### Length and weight parameters follow the one you gave last

A size can be given either as a weight (`w_mat`, `w_max`, …) or as the length it
converts to (`l_mat`, `l_max`, …). Mizer used to derive the weight from the
length whenever both were present, so on a model specified by lengths a weight
could not be set at all: the value you assigned was replaced on the spot by the
one calculated from the unchanged length.

Both now follow one rule: **the one you gave last wins, and if you gave both at
the same time the weight wins.** The other is set to match, so the two never
disagree, and mizer warns, naming the species, when it changes a length to match
a weight it disagrees with.

```r
params <- newMultispeciesParams(sp)   # sp specifies l_mat, a and b

# Used to be silently undone, now it takes effect and l_mat follows
species_params(params)$w_mat[1] <- 100
```

The rule is applied when a data frame is assigned into a model, which is when
mizer can tell which values changed. A data frame you have taken out of a model
and are editing on its own is left exactly as you write it — the conversions,
checks and warnings happen on assignment. One that was never in a model, for
example one passed to `validSpeciesParams()`, carries no such history, so a
length and a weight that disagree there count as given at the same time and the
weight wins.

`species_params<-()` and `given_species_params<-()` apply the rule identically
(#490). Previously only `species_params<-()` did, so the same edit made through
`given_species_params<-()` was discarded — a maturity weight differing by up to
73% — and the given species parameters were left permanently inconsistent, which
made mizer repeat the "not consistent" warning at every later parameter change.
If you worked around this by setting the weight and the length together, you can
now set either one on its own.

### A change that cannot take effect now warns

If you have set a rate array by hand — `metab(params) <- ...`, or any setter
called with an explicit array — mizer no longer calculates that rate from the
species parameters. Changing one of the species parameters that used to feed it
therefore changes the species parameter table but not the model. Mizer now
warns when that happens:

```r
params <- NS_params
metab(params) <- metab(params)      # freeze the metabolic rate

sp <- species_params(params)
sp$ks <- sp$ks * 2
species_params(params) <- sp
#> Warning: Your change to the species parameter `ks` has not taken effect
#> because the metabolic rate has been set manually and so is no longer
#> calculated from the species parameters. Call
#> `setMetabolicRate(params, reset = TRUE)` if you want the metabolic rate to
#> be calculated from the species parameters again.
```

Previously this was silent: the rate setter did emit a message, but
`species_params<-()` and `given_species_params<-()` run `suppressMessages()`
over the recalculation to quieten the routine chatter, so the message never
reached you (#489).

**How this affects existing code:** nothing about the model changes — the
change did not take effect before either, you just were not told. But code that
runs under `options(warn = 2)`, or a test using `expect_message()` where the
report now arrives as a warning, will need adjusting. To act on the warning,
either put the rate back under the control of the species parameters with the
`reset = TRUE` call it names, or set the rate array itself instead of the
species parameter. To silence it, set `options(mizer_info_level = 0)`, which
also covers the functions that take no `info_level` argument, or pass
`info_level = 0` to the one call you want quiet.

The resource works the same way. `resource_params(params)$kappa <- ...` on a
model whose `resource_capacity()` you had set by hand used to change the stored
`kappa` and nothing else, without saying anything at all; it now warns.

`species_params<-()` also now gives the warning that `given_species_params<-()`
already gave when a change is ignored because another parameter takes
precedence over it — a change to `f0` when `gamma` has been given, to `fc` when
`ks` has, or to `age_mat` when `h` has.

This is the counterpart of the 3.2 change described under *Frozen arrays are
protected from incidental balancing* below: there, mizer keeps a frozen array
instead of overwriting it; here, it tells you when a parameter change is
ignored because an array is frozen.

### One report, one switch

Nearly everything mizer says while building or changing a model now goes
through the same mechanism, including the reports in `steady()`,
`projectToSteady()`, `matchYields()`, `validParams()`, `setInteraction()`,
`setReproduction()`, `setResource()`, `newTraitParams()`,
`newSingleSpeciesParams()` and `plotYieldObservedVsModel()`. Two consequences
for existing code:

- **`info_level = 0` now means silence.** Reports that were plain `message()`
  calls ignored `info_level` altogether and appeared anyway; they no longer do.
  If your code relied on seeing one of them, drop the `info_level = 0`.
- **Reports are collected and given at the end of the call**, one message and
  one warning rather than a stream. A test doing `expect_message()` on an
  individual report inside a longer call may need adjusting, and the text now
  arrives with any others in the same message.

### `$` on a parameter table no longer partially matches

`$` on a `species_params` or `gear_params` table now matches column names
exactly. Partial matching was silently returning the wrong parameter: in a model
without the length-weight parameters `a` and `b`,

```r
species_params(NS_params)$a   # used to return the `alpha` column
species_params(NS_params)$b   # used to return the `beta` column
```

complete with per-species names, so code converting weights to lengths got the
assimilation efficiency and the preferred predator/prey mass ratio instead.
Writing was never partially matched (`sp$b <- 3` always created a new column
`b`), so reads and writes disagreed about what `$b` meant.

A name that is not a column now gives `NULL`. If the name would have partially
matched a single column, you also get a warning naming that column. So
`is.null(species_params(params)$foo)` is now a reliable test for whether a
parameter is present, and any code that was relying on the abbreviation should
spell the column out in full (#487).

### Quadrature fixes under `second_order_w`

Two diagnostics were applying the prey-bin quadrature twice when second-order
bin-averaging was switched on with `second_order_w()`, and are now consistent
with `getEncounter()`:

- **`getDiet(proportion = FALSE)`** was uniformly too large by a factor
  `(1 + beta) / 2`, where `beta` is the grid ratio — 9.7% for `NS_params`.
  Summing the diet over prey now reproduces
  `getEncounter() * (1 - getFeedingLevel())` under both schemes (#474).
  `getDiet(proportion = TRUE)`, the default, was unaffected: the factor was
  uniform and divided out.
- **`getTrophicLevel()`** built its numerator and denominator from different
  quadratures, so reported trophic levels were off by up to 0.06. A predator
  whose prey all have trophic level 1 now comes out at exactly 2 under both
  schemes (#474).

Models on the default (first-order) scheme are unchanged. If you have published
absolute diet values or trophic levels computed under `second_order_w`, they need
recomputing.

### `getProportionOfLargeFish()` on a `MizerParams` object was wrong

The `MizerParams` method multiplied the species x size abundance array by the
vector of weights, which R recycles down the columns of the array rather than
along the size axis, so every species but the first was weighted by the wrong
sizes. Only the `MizerParams` method was affected; the `MizerSim` method was
always right, and the two now agree when applied to the same state (#494). Any
Large Fish Index computed from a `MizerParams` object in a model with more than
one species needs recomputing.

### `getN()` respects the quadrature scheme at the ends of the size range

`getN(params, min_w = ...)` now bin-averages the size-range window when
second-order bin-averaging is switched on with `second_order_w()`, so the bin
straddling `min_w` or `max_w` contributes only partially — as `getBiomass()`
already did. Numbers over a restricted size range therefore change slightly
under `second_order_w`; over the full size range, and on the default
first-order scheme, nothing changes (#494).

### `w_min` survives a rebuild of the species parameters

`w_min` is now part of `given_species_params`, so the `min_w` argument to
`newMultispeciesParams()` and `emptyParams()` is preserved across any operation
that rebuilds the species parameters. Previously a `given_species_params<-`
round-trip silently reset `w_min` to 0.001 when `min_w` was smaller, and emitted
a spurious warning when it was larger (#460). Code that set a small `min_w` and
worked around the reset — or that unknowingly ran on the reset grid — now gets
the size grid it asked for, and results change accordingly.

### `compareParams()` compares small parameters properly

`compareParams()` now uses a relative tolerance for species parameters, so
small-magnitude parameters such as `gamma` (~1e-8) are no longer treated as equal
when they differ by up to ~10%. Comparisons that previously reported two models
as identical may now report differences — those differences were always there.

### Defaults for `gamma` and `f0` ignore a hand-set search volume

`get_gamma_default()` works out how much energy is available to a predator by
giving it a search volume coefficient of 1. It used to obtain that search volume
by calling `setSearchVolume()`, which refuses to recalculate a `search_vol`
array you have set by hand — so mizer's own internal call was blocked along with
yours, and the available energy was measured with *your* array. The resulting
`gamma` was wrong by whatever factor separated your array from the unit-gamma
one, which in a realistic model is many orders of magnitude. `get_f0_default()`,
the inverse, had the same problem. Both now build the search volume they need
directly from the species parameters (#488).

```r
sv <- search_vol(params)
search_vol(params) <- sv * 10          # freeze the search volume

sp <- species_params(params)
sp$gamma <- NA
species_params(params) <- sp           # ask mizer to recalculate gamma

species_params(params)$gamma
#> Used to come back ~1e9 times too large; now the same value you would
#> get without the frozen search volume.
```

If you have a model in which you set `search_vol` by hand and then let mizer
fill in a missing `gamma` or `f0`, that model's species parameters were wrong
and change with this release. Note that the recalculated `gamma` still has no
effect on the model while the search volume stays frozen — mizer now warns you
about that separately, see "A change that cannot take effect now warns" above.
Call `setSearchVolume(params, reset = TRUE)` to put the search volume back under
the control of the species parameters.

### `setParams()` rejects arguments it does not use

`setParams()` passes its `...` on to the rate setters, each of which declares
its own `...` as unused. Any argument that none of them recognises was
therefore accepted and ignored without a word. It is now an error, and the
error says where the argument belongs when it belongs somewhere:

```r
setParams(params, metabolic = 99)        # was: silently ignored
setParams(params, resource_rate = 5)     # was: silently ignored
```

The resource case is the one most likely to have bitten: `setParams()` never
called `setResource()`, so no resource argument ever reached the model, and the
deprecation warnings for `setResource(r_pp)` and `setResource(kappa)` used to
recommend `setParams(resource_rate)` and `setParams(resource_capacity)`, which
do nothing. Use `setResource()` for all of these:

```r
params <- setResource(params, resource_rate = 5)
```

Likewise `gear_params` goes to `gear_params<-()`, and `second_order_w` and
`use_predation_diffusion` to their own assignment functions. `setResource()`
itself now applies the same check to its own `...`, so a misspelled resource
argument errors there too.

If your code errors here, the argument was having no effect before, so removing
it changes nothing; moving it to the right function changes the model, and that
is the change you had intended all along.

Two related tidy-ups: `reset` is now a documented argument of `setParams()`
(it was already forwarded through `...`, undocumented) and still thaws every
rate array that `setParams()` sets; and `setExtDiffusion()` is now listed
among the setters that `setParams()` calls, which it always did.

<!-- agent-only -->

Diagnostic: if a user reports that a resource change "did not take", check
whether they went through `setParams()`. Before this release the call was
accepted, so there is no error in their logs and the model simply kept its old
resource. `resource_rate(params)` before and after their call is the quickest
confirmation.

<!-- /agent-only -->

### One name for each stored rate array

Eleven accessors that read a rate array back out of a `MizerParams` object had
two interchangeable names. The bare name is now the one to use — it is the one
that also has a replacement function, so the pair reads the same way in both
directions (`catchability(params)` and `catchability(params) <- value`). The
`get`-prefixed names are soft-deprecated and warn:

| Deprecated | Use instead |
|---|---|
| `getCatchability()` | `catchability()` |
| `getSelectivity()` | `selectivity()` |
| `getInitialEffort()` | `initial_effort()` |
| `getPredKernel()` | `pred_kernel()` |
| `getSearchVolume()` | `search_vol()` |
| `getMaxIntakeRate()` | `intake_max()` |
| `getMetabolicRate()` | `metab()` |
| `getExtMort()` | `ext_mort()` |
| `getExtEncounter()` | `ext_encounter()` |
| `getMaturityProportion()` | `maturity()` |
| `getReproductionProportion()` | `repro_prop()` |

Nothing breaks: the old names still return exactly the same value. Renaming is
a search and replace.

The `get` prefix now means one thing — a function that *calculates* something
from the current state of a model, like `getEncounter()`, `getFMort()` or
`getBiomass()`. The functions above only hand back a value that is already
stored in the object.

<!-- agent-only -->

The `get` forms are no longer S3 generics; they are plain functions that warn
and forward. Dispatch still works for a custom class, but on the bare name, so
an extension that defined `getExtMort.MyClass` must rename its method to
`ext_mort.MyClass`. A method on the bare name has always worked and keeps
working through both names.

When a user's `bin_average` diagnostic disagrees with the rate functions, check
whether they reached for `pred_kernel()` (formerly `getPredKernel()`) rather
than `encounter_kernel()` — the rename does not change that distinction, but it
makes the two names look more alike than they used to.

<!-- /agent-only -->

## Upgrading from mizer 3.1 to 3.2

### `species_params<-()` now detects and protects changes

Previously, modifying species parameters via `species_params<-()` updated the values in the model but bypassed `given_species_params()`. This meant that your changes were not protected, and any subsequent recalculation of defaults (for example, by a call to `given_species_params<-()`) would overwrite your custom values. Furthermore, changing a parameter like `w_inf` via `species_params<-()` did not automatically trigger a recalculation of downstream parameters like `w_mat` or `w_max`.

Now, `species_params<-()` intelligently diffs the new data frame against the old one to detect exactly which parameters you have changed. It automatically records those changed parameters in `given_species_params`, protecting them from future overwrites, and immediately recalculates any downstream defaults based on your changes. 

**How this affects existing code:**

1. If your existing code used `species_params<-()` to update a core parameter like `w_inf` and you expected `w_mat` or `w_max` to remain frozen at their old values, you will now see them automatically recalculate. If you wish to freeze downstream parameters, you must provide their frozen values explicitly in the same update.

2. If your code computes custom parameters and saves them via `species_params<-()`, those parameters will now be preserved and survive future recalculations.


### Setting resource parameters

Two related changes affect how you modify the resource size spectrum. Together
they make the resource scalars behave like the species parameters: a scalar is
an input, and the size-dependent arrays are computed from it.

#### Assigning to `resource_params()` now updates the resource arrays

Previously, assigning to `resource_params()` — or to one of its components, such
as `resource_params(params)$kappa <- ...` — only stored the new scalar values.
The size-dependent carrying capacity (`cc_pp`) and replenishment rate (`rr_pp`)
were left unchanged until you next called `setResource()`.

Now these assignments immediately rebuild the arrays from the scalars, exactly as
`species_params()<-` rebuilds the species rates:

- `kappa`, `lambda` and `w_pp_cutoff` rebuild the carrying capacity;
- `r_pp` and `n` rebuild the replenishment rate.

Arrays that you have set by hand are left untouched (see *Frozen arrays* below).

If your code changed a resource scalar and then called `setResource()` to apply
it, nothing breaks — you can drop the now-redundant `setResource()` call. If you
changed a resource scalar and relied on the arrays *not* changing until later,
review that code.

#### Assigning to `resource_params()` does not balance the resource

*Balancing* means adjusting the rate and capacity together so that the resource
replenishes at exactly the rate at which it is consumed, keeping it at its steady
state. Assigning to `resource_params()` rebuilds the arrays from the scalars but
does **not** balance, so the resource steady state generally shifts.

Balancing is now solely a feature of `setResource()`. To change a resource
coefficient *and* keep the resource balanced, call `setResource()` rather than
assigning to `resource_params()`:

```r
# Rebuild the capacity from a new coefficient and rebalance the rate,
# so the steady state is preserved:
params <- setResource(params, resource_capacity = new_kappa)

# Likewise, set a new rate coefficient and rebalance the capacity:
params <- setResource(params, resource_rate = new_r_pp)
```

#### The resource setters gained a `balance` argument

`resource_rate<-`, `resource_capacity<-`, `resource_level<-` and
`resource_dynamics<-` still balance by default (unchanged behaviour), but they
now accept a `balance` argument so you can switch balancing off:

```r
# Set the capacity but leave the rate untouched (do not rebalance):
resource_capacity(params, balance = FALSE) <- my_capacity
```

#### Frozen arrays are protected from incidental balancing

When you set the size dependence of the resource capacity or the resource rate
by hand (by assigning a full vector rather than a scalar), mizer marks it "set
manually" — it is *frozen* and will not be recomputed from the resource
parameters. Previously, an operation that re-balanced the resource *without*
being given a replacement rate or capacity — for example changing only
`resource_dynamics`, or calling `setResource()` with neither a rate nor a
capacity — would silently overwrite such a frozen array. It is now kept, and a
warning is issued instead.

To deliberately recompute a frozen array from the resource parameters, pass
`reset = TRUE` to `setResource()`.

### The `species_params` data frame is now an S3 subclass

The `species_params` data frame now has class `c("species_params",
"data.frame")` (and `gear_params` similarly). It behaves like an ordinary data
frame, but subsetting and subassignment go through class-preserving S3 methods
and can trigger reactive re-validation and conversions (for example filling in a
weight from a length). Code that relied on `class(species_params(params))` being
exactly `"data.frame"`, or that stripped attributes with the assumption of a
plain data frame, may need adjusting. When you need a plain frame, coerce
explicitly with `as.data.frame()`.

### Accessing a column with `$` now returns a named vector

Extracting a single column from a `species_params` or `gear_params` object with
`$` now returns a vector named by species (or by `"species, gear"` for
`gear_params`):

```r
species_params(params)$w_mat
#>   Sprat  Herring      Cod
#>    ...      ...      ...
```

The values are unchanged, but the names are new. This is convenient for
identifying entries, but code that compared such a vector with `identical()` to
an unnamed vector, or that used it as-is where names matter (for example as
row/column names elsewhere), may behave differently. Strip the names with
`unname()` if you need the old behaviour. The `species` column itself is
returned unnamed.

### Setting `sel_func` adds the required argument columns

Assigning a selectivity function name to a `gear_params` object now
automatically adds the argument columns that the function needs (as `NA`),
ready to be filled in:

```r
gp$sel_func <- "sigmoid_length"
# gp now has l25 and l50 columns, both NA
```

Previously these columns had to be added by hand. Code that checks which columns
are present in `gear_params`, or that expected setting `sel_func` to leave the
column set unchanged, will now see the extra columns (#431).

### Passing a data frame to `species_params()` / `given_species_params()` now validates it

Calling `species_params()` or `given_species_params()` on a plain data frame now
runs the same validation and defaults that `validSpeciesParams()` and
`validGivenSpeciesParams()` apply, rather than only checking for misspellings and
converting lengths to weights. `species_params(df)` fills in the default columns
(`w_max`, `alpha`, `n`, `p`, `interaction_resource`, `z_ext`, and the rest), and
`given_species_params(df)` applies the consistency corrections (for example
clamping `w_mat` below `w_inf`), derives `w_inf` from `w_max`/`w_repro_max` when
it is absent, and now stops if the frame has duplicate species rows. Models built
or modified through `newMultispeciesParams()`, `setParams()` and the
`species_params()<-` / `given_species_params()<-` setters are unaffected, because
those already ran this validation. Only code that called the two accessors
directly on a bare data frame will see the extra columns and stricter checks
(#432).

### Printing of mizer array objects shows the values

`print()` on the array objects returned by the rate getters (`ArraySpeciesBySize`,
`ArrayTimeBySpecies`, `ArrayResourceBySize`, `ArrayTimeByResourceBySize` and
`ArrayTimeBySpeciesBySize`, as returned by `getEncounter()`, `getBiomass()`,
`getFMort()`, `NResource()` and similar) now truncates the output instead of
flooding the console with all the array entries. If your code or reports relied
on the old printed format, use `as.data.frame()` to go back to the full output.

### Upper boundary condition at `w_max`

The size-spectrum solver now holds the abundance at zero above each species'
maximum size `w_max`. Without diffusion this happens automatically and results
are unchanged. With diffusion switched on this change stops a small amount of
density leaking to sizes above `w_max`, so results there change slightly. See
`vignette("numerical_details")`.

### Extension packages: dynamic marker classes

If you develop a mizer extension, an installed extension package is now
recognised as a dispatching extension from the S3 methods it registers for its
marker class (for example `getEncounter.mizerMR`), rather than only from a
statically defined S4 marker class. You can now omit the static
`setClass("mizerFoo", contains = "MizerParams")` and let mizer create the marker
class dynamically. This lets two independently developed extensions be chained
in either load order. See `vignette("creating-extension-packages")`.

## Upgrading from mizer 3.0 to 3.1

Version 3.1 leaves default results unchanged from 3.0 unless you opt in to the
new experimental second-order-in-size scheme. The changes below can still affect
existing code in specific situations.

### Maximum-size species parameters clarified

The maximum-size parameters have been given clearer, separate roles (#325):

- `w_inf`, the von Bertalanffy asymptotic size, is now the primary maximum-size
  parameter and is used as the default for `w_repro_max` (the size at which a
  mature individual invests all its energy in reproduction) and for `w_mat`.
- `w_max` is now purely a computational boundary — it sets the size grid and the
  plot range — and defaults to `1.5 * w_inf`.
- The default external mortality parameter `z0` is now computed from `w_inf`
  rather than `w_max`, so the computational boundary `w_max` no longer feeds into
  any model parameter.

Existing models and scripts are unaffected: if `w_inf` is not supplied it is
taken from `w_repro_max` or `w_max`, so old objects behave as before. However,
**new models built from the defaults may differ from 3.0.0**. If you build models
from scratch, check that `w_inf`, `w_max` and `w_repro_max` mean what you intend.

### `getTrophicLevel()` gives the resource a size-dependent trophic level

`getTrophicLevel()` and `getTrophicLevelBySpecies()` now assign the resource a
size-dependent trophic level,
$T_R(w) = \max(1,\, 1 + \log(w / w_R) / \log(\beta_R))$, instead of treating the
resource as trophic level 0. The new `w_R` and `beta_R` arguments control this.
Trophic levels computed with these functions will therefore be higher than
before. Set the arguments explicitly if you need to reproduce old numbers.

### Bug fixes that change results

Several fixes correct earlier behaviour and so change output:

- **`summary()` of a `MizerSim`** now reports the fishing effort that was used
  during the simulation, rather than the model's `initial_effort`. Gears whose
  effort varied over time show the mean, flagged with a note giving the range.
  The printed summary therefore differs for simulations run with time-varying
  effort.
- **`MizerSim` method for `plotDiet()`** introduced in version 3.0 simply
  plotted the diet at the initial time of the simulation. Now `plotDiet()` for
  a `MizerSim` accepts a `time_range` argument. The diet is now computed from
  the *simulated* abundances at the requested times, defaulting to the final
  saved step, rather than the initial one (#357).
- **Other components and `t_save`.** `project()` was advancing the abundances of
  other components (set via `setComponent()`) only once per *saved* time step
  instead of once per `dt` step. They are now integrated with the same `dt` as
  the consumer and resource spectra, so results for models with other components
  no longer depend on `t_save`.
- **Time-varying effort in `getRDI()`, `getRDD()`, `getFlux()`.** On a
  `MizerSim` object these now use the simulated time-varying effort rather than
  the initial effort, so they change for simulations with varying effort (#370).
- **`plotCDF()` / `plotlyCDF()` bin placement.** Each cumulative value is now
  plotted at its bin's *upper* edge, correcting a one-bin offset.
  The curves shift by one bin compared with 3.0 (#383).
- **`distanceMaxRelRDI()`.** Now returns `Inf` instead of `NaN` when a previous
  RDI is zero, so `projectToSteady()` no longer mistakes a `NaN` distance for
  convergence. Convergence behaviour can therefore differ in edge cases.

### Second-order methods advance the resource at the midpoint

If you use `project()` with `method = "predictor_corrector"` (or the new
`method = "tr_bdf2"`), the resource and the other components are now advanced
with midpoint rates rather than the start-of-step value, so that they reach the
same second-order accuracy in time as the consumer spectra. Results from these
methods therefore differ slightly from 3.0. The default `method = "euler"` and
the steady states are unchanged.

### Opting in to the second-order-in-size scheme

3.1 adds an optional, experimental second-order-accurate finite-volume scheme in
the size variable, controlled by the new `second_order_w` slot. It is **off by
default**, so default results are byte-identical to 3.0. If you switch it on
(via `second_order_w()<-` or the `second_order_w` argument of the `new...Params()`
constructors), size-integrated diagnostics and the resource spectrum shift by
$O(\Delta w)$, so a calibrated model may need recalibrating. See `?second_order_w`
and the "Numerical Details" vignette.

## Upgrading from mizer 2.5.4 to 3.0

Version 3.0 is a large release. Most new capabilities are additive and off by
default, but there are several renamed arguments, deprecations and behavioural
changes that can affect existing code.

### Renamed arguments and changed defaults (breaking changes)

- **First argument of `plotBiomass()`, `plotYield()`, `plotYieldGear()`** (and
  their `MizerSim` methods and `plotly*` wrappers) is renamed from `sim` to
  `object`, for consistency with the other plot generics. Calls that passed the
  simulation by name, `plotBiomass(sim = my_sim)`, must become
  `plotBiomass(object = my_sim)`. Positional calls are unaffected.
- **`plotBiomassObservedVsModel()` / `plotlyBiomassObservedVsModel()`** now
  default to `ratio = FALSE` for all object types. Calls that relied on the
  previous ratio plot must set `ratio = TRUE` explicitly.
- **`plotDiet()` no longer accepts a `time_range` argument.** Remove it from your
  calls. (In 3.1 a `time_range` argument returns for the `MizerSim` method — see
  above.)
- **Dimnames of `getMort()` and `getPredRate()`** arrays are now `sp` and `w`
  (matching `getFMort()` and the other rate getters). Code that referred to the
  old dimnames by name must be updated.

### Rate getters return classed array objects

Functions that return arrays of the form (species × size), (time × species) or
(time × species × size) now attach extra attributes and an S3 class
(`ArraySpeciesBySize`, `ArrayTimeBySpecies` or `ArrayTimeBySpeciesBySize`). The
numeric values and ordinary matrix behaviour (arithmetic, subsetting) are
unchanged, but the extra class and attributes mean that a strict comparison such
as `identical(getMort(params), old_value)` can now report a difference where the
numbers agree. Use `unclass()`, or compare with `all.equal()` on the values, if
you need to ignore the class. These objects also carry `print()`, `summary()`,
`plot()` and `as.data.frame()` methods, so printing them looks different from a
bare matrix.

### `setInitialValues()` is deprecated

`setInitialValues()` is deprecated. Replace

```r
params <- setInitialValues(params, sim)
```

with

```r
params <- finalParams(sim)
```

or, when averaging over a time range, with
`getParams(sim, time_range, geometric_mean)`. This reflects a shift in
interpretation: a `MizerParams` object now represents not just the model
specification but also its current state (the abundances), which can be
extracted from a simulation with `getParams()`, `finalParams()` and
`initialParams()`.

### Growth can no longer be negative

Growth is now forced to be non-negative, preventing unphysical shrinkage. In any
model where the energy available for growth used to go negative (for example a
strongly food-limited large individual), growth is now clamped at zero instead,
so projected size spectra can differ from 2.5.4. No warning is issued when growth
stops at or after the maturity size.

### `project()` timing and effort handling

- **Inherited `dt` and `method`.** When `project()` is called on an existing
  `MizerSim` object, `dt` and `method` now default to the values stored in the
  simulation's new `sim_params` slot. If you pass values that differ from the
  stored ones, a warning is issued. To use different settings deliberately, pass
  them explicitly and expect the warning.
- **`t_max` / `t_save` with an effort array.** These arguments are now respected
  even when an effort array is supplied (#231). With `t_max` the simulation
  extends beyond the times in the effort array using the last known effort; with
  `t_save` the save frequency is controlled independently, interpolating effort
  as needed. Simulations that previously derived their length or save times
  solely from the effort array may now produce a different set of saved steps.
- **State at `t_max` always saved.** `project()` now warns when `t_max` is not a
  multiple of `t_save` and ensures the state at `t_max` is saved even if the
  final interval is shorter than `t_save` (#341). The returned simulation may
  therefore contain one extra saved time step compared with 3.0.

### `plot()` and `summary()` are now S3 methods

The `plot()` and `summary()` methods for `MizerParams`, `MizerSim` and the mizer
array classes are now registered as S3 methods rather than S4 methods, so
`plot()` and `summary()` stay plain S3 generics when mizer is loaded. This avoids
interfering with S4 dispatch in other packages, but code that relied on
`plot`/`summary` being S4 generics (for example via `selectMethod()` or
`getMethod()`) needs adjusting.

### Bug fixes that change results

- **`getMeanMaxWeight()`** now applies the species selector to the denominator as
  well, so its values change when a subset of species is selected.
- **`plotSpectra()` axis limits.** It no longer forces the y-axis lower limit to
  `1e-20` (it auto-scales to the data) and, when `resource = FALSE`, it uses
  `min(params@w)` rather than `min(params@w) / 100` as the default lower size
  limit. Plots therefore look different.
- **`getFMort()` on a `MizerSim`** was silently dropping the component names from
  `n_other`, breaking rate functions that access `n_other` by name (e.g.
  `n_other[["resource"]]`); it now preserves them.
- **`getFMort.MizerSim()`** now passes the time argument `t` to user-defined
  fishing-mortality functions, so a time-dependent fishing function now sees the
  correct time.

### Predation diffusion is available but off by default

3.0 adds a diffusion term to the growth dynamics, controlled by the new
`use_predation_diffusion` slot. It defaults to `FALSE`, preserving the behaviour
of earlier mizer, so existing models are unchanged unless you switch it on with
`use_predation_diffusion(params) <- TRUE`. Likewise the new species parameters
`z_ext`, `d`, `E_ext` and `D_ext` for external mortality, encounter and diffusion
all default to values that leave the model unchanged.
